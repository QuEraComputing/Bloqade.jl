
"""
    struct ACFETEvolution
        ACFETEvolution(reg::AbstractRegister, start_clock, end_clock, h, algo; step_size = 1e-7, kw...)

Create a `ACFETEvolution` object that describes a time evolution
using commutation-free Magnus method with s exponentials with adaptive size
[1] https://arxiv.org/pdf/1102.5071.pdf (eq.58-61) 
[2] SAIM: M2AN 53 (2019) 197–218 A posteriori error estimation for Magnus-type integrators W. Auzinger et.al.

# Arguments

- `reg`: a register, should be a subtype of `AbstractRegister`.
- `start_clock`: the start clock of time evolution.
- `end_clock`: the end clock of time evolution. 
- `h`: a hamiltonian expression.
- `algo`: the algorithm (different orders to use), should be a subtype of `CFETTables`. default is `CFET2_1`

# Keyword Arguments

- `step_size`: intial step size, default is `1e-7`.
- `progress`: show progress bar, default is `false`.
- `progress_name`: progress bar name, default is `"emulating"`.
- `normalize_step`: normalize the state every `normalize_step`.
- `normalize_finally`: wether normalize the state in the end of evolution, default is `true`.
- `tol`: tolerance of the exponential time propogator evaluation method, default is `1e-7`
- `expmv_backend`: the backend for evaluate exponential time propogator, default is `expmv!` (other option: `expm_multiply!`).

# Examples

The following is the simplest way of using `ACFETEvolution`
via [`emulate!`](@ref). For more advanced usage, please refer
to documentation page [Emulation](@ref emulation).

```jldoctest
julia> using Bloqade

julia> r = zero_state(5)
ArrayReg{2, ComplexF64, Array...}
    active qudits: 5/5
    nlevel: 2

julia> atoms = [(i, ) for i in 1:5]
5-element Vector{Tuple{Int64}}:
 (1,)
 (2,)
 (3,)
 (4,)
 (5,)

julia> h = rydberg_h(atoms; Ω=sin)
nqubits: 5
+
├─ [+] ∑ 5.42e6/|x_i-x_j|^6 n_i n_j
└─ [+] Ω(t) ⋅ ∑ σ^x_i


julia> prob = ACFETEvolution(r, 0.0, 0.1, h, CFET2_1());

julia> emulate!(prob); # run the emulation
```
"""
mutable struct ACFETEvolution{Reg<:AbstractRegister,T<:Real,H<:Hamiltonian,ALGTBL<:CFETTables} <: ADEvolver
    reg::Reg
    start_clock::T
    end_clock::T
    step_size::T # this store the current step_size
    hamiltonian::H
    options::KrylovOptions
    alg_table::ALGTBL

    function ACFETEvolution{Reg,T,H,ALGTBL}(reg, start_clock, end_clock, step_size,hamiltonian, options, algo) where {Reg,T,H,ALGTBL}
        start_clock ≥ 0 || throw(ArgumentError("start clock must not be negative"))
        end_clock ≥ 0 || throw(ArgumentError("end clock must not be negative"))

        end_clock ≥ start_clock || throw(ArgumentError("end clock must be larger than start clock"))
        step_size ≤ end_clock-start_clock || throw(ArgumentError("initial step size cannot be larger than the whole evolution time"))
        return new{Reg,T,H,ALGTBL}(reg, start_clock, end_clock, step_size, hamiltonian, options, algo)
    end
end

"""
    ACFETEvolution(reg, start_clock, end_clock, hamiltonian, options)

Create a `ACFETEvolution` object.

# Arguments

- `reg`: a register object.
- `start_clock`: start clock of the evolution.
- `end_clock`: end clock of the evolution.
- `hamiltonian`: low-level hamiltonian object of type [`Hamiltonian`](@ref).
- `options`: options of the evolution in type [`KrylovOptions`](@ref).
"""
function ACFETEvolution(reg, start_clock, end_clock, step_size,hamiltonian, options, algo)
    return ACFETEvolution{typeof(reg),typeof(start_clock),typeof(hamiltonian),typeof(algo)}(
        reg,
        start_clock,
        end_clock,
        step_size,
        hamiltonian,
        options,
        algo,
    )
end

function Adapt.adapt_structure(to, x::ACFETEvolution)
    return ACFETEvolution(adapt(to, x.reg), x.start_clock, x.end_clock, x.step_size, adapt(to, x.hamiltonian), x.options, x.alg_table)
end

function ACFETEvolution(reg::AbstractRegister, start_clock, end_clock , h, algo::CFETTables = CFET2_1(); step_size=1e-7, kw...)
    all(≥(0), clocks) || throw(ArgumentError("clocks must not be negative"))
    options = from_kwargs(KrylovOptions; kw...)
    P = real(eltype(statevec(reg)))
    T = isreal(h) ? P : Complex{P}
    
    return ACFETEvolution(reg, start_clock, end_clock, step_size, Hamiltonian(T, h, space(reg)), options, algo)
end

## given Hamiltonian, current time `t` and duration `dt`, plus a CFETTable, 
# output generator Ω at certain ETstep (exponential time propogator) 
## t + xs[i]dt
function __construct_Ω(h::Hamiltonian, t::Real, dt::Real, Tbl::CFETTables, ETStep::Int)
    
    gs = Tbl.Gs[ETStep]
    xs = Tbl.xs 
    
    fs = gs[1]*h(t + xs[1]*dt).fvals

    for i in 2:length(gs)      
        fs += gs[i]*h(t + xs[i]*dt).fvals
    end

    return SumOfLinop{LinearAlgebra.Hermitian}(fs, h.ts)
end

function __construct_dΩ(h::Hamiltonian, t::Real, dt::Real, Tbl::CFETTables, ETStep::Int)
    
    gs = Tbl.Gs[ETStep]
    xs = Tbl.xs 
    
    fs = gs[1]*xs[1]*derivative.(h.fs,t + xs[1]*dt)

    for i in 2:length(gs)      
        fs += gs[i]*xs[i]*derivative.(h.fs,t + xs[i]*dt)
    end

    return SumOfLinop{LinearAlgebra.Hermitian}(fs, h.ts)
end

@inline function scale_factor(αmin, αmax, α, ϵ, tol, p)
    return min(αmax, max(αmin, α*(tol/ϵ)^(1/(p+1))))
end

## here, since generic algo for defect is unclear, we specialize it for each support types:
function emulate_step!(prob::ACFETEvolution{<:Any,<:Any,<:Any,ALGTBL<:CFET2_1}, step::Int, clock::Real, tol::Real)
    state = statevec(prob.reg)
    Ham = prob.hamiltonian

    duration = prob.step_size # get duration 

    #construct Ωi:
    Ω1 = -im*__construct_Ω(Ham, clock, duration, prob.alg_table, 1)
    dΩ1 = -im*__construct_dΩ(Ham, clock, duration, prob.alg_table, 1)


    # perform evolution:
    prob.options.expmv_backend(duration, Ω1, state; prob.options.tol)


    ## calculate defect operator
    Bv = similar(state)
    dBv = similar(state)
    tmp = similar(state)

    mul!(Bv, Ω1, state)
    mul!(dBv, dΩ1, state)
    dest = Bv .+ duration*dBv

    mul!(tmp, Ω1, dBv)
    mul!(dBv, dΩ1, Bv)
    Bv = tmp

    ## [B,B']v
    dest .+= (duration^2/2)*(Bv - dBv)

    ## -A(t+τ)v
    mul!(tmp, Ham(clock+duration), state)
    dest .-= tmp

    ϵ = norm(dest)

    ## calculate and update next step size:
    prob.step_size = duration*scale_factor(0.25, 4.0 , 0.9, ϵ, tol, p)


    # do we need this normalization? 
    if mod(step, prob.options.normalize_step) == 0
        normalize!(prob.reg)
    end

    if prob.options.normalize_finally && step == length(prob.durations)
        normalize!(prob.reg)
    end
    
    return prob
end

function emulate_step!(prob::ACFETEvolution, step::Int, clock::Real, duration::Real)
    error("unsupported CFET algo for adaptive step-size")
end

