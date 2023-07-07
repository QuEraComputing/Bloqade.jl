abstract type CFETTables end

"""
    struct CFETEvolution
        CFETEvolution(reg::AbstractRegister, clocks, h, algo; kw...)

Create a `CFETEvolution` object that describes a time evolution
using commutation-free Magnus method with s exponentials.
https://arxiv.org/pdf/1102.5071.pdf (eq.58-61) 


# Arguments

- `reg`: a register, should be a subtype of `AbstractRegister`.
- `clocks`: the clocks of this time evolution at each step.
- `h`: a hamiltonian expression.
- `algo`: the algorithm (different orders to use), should be a subtype of `CFETTables`. default is `CFET4_2`

# Keyword Arguments

- `progress`: show progress bar, default is `false`.
- `progress_name`: progress bar name, default is `"emulating"`.
- `normalize_step`: normalize the state every `normalize_step`.
- `normalize_finally`: wether normalize the state in the end of evolution, default is `true`.
- `tol`: tolerance of the exponential time propogator evaluation method, default is `1e-7`
- `expmv_backend`: the backend for evaluate exponential time propogator, default is `expmv!` (other option: `expm_multiply!`).

# Examples

The following is the simplest way of using `CFETEvolution`
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


julia> prob = CFETEvolution(r, 0.0:1e-2:0.1, h, CFET4_2());

julia> emulate!(prob); # run the emulation
```
"""
struct CFETEvolution{Reg<:AbstractRegister,T<:Real,H<:Hamiltonian,ALGTBL<:CFETTables} <: Evolver
    reg::Reg
    start_clock::T
    durations::Vector{T}
    hamiltonian::H
    options::KrylovOptions
    alg_table::ALGTBL

    function CFETEvolution{Reg,T,H,ALGTBL}(reg, start_clock, durations, hamiltonian, options, algo) where {Reg,T,H,ALGTBL}
        start_clock ≥ 0 || throw(ArgumentError("start clock must not be negative"))
        all(≥(0), durations) || throw(ArgumentError("durations must not be negative"))
        return new{Reg,T,H,ALGTBL}(reg, start_clock, durations, hamiltonian, options, algo)
    end
end

"""
    CFETEvolution(reg, start_clock, durations, hamiltonian, options)

Create a `CFETEvolution` object.

# Arguments

- `reg`: a register object.
- `start_clock`: start clock of the evolution.
- `durations`: list of durations at each time step.
- `hamiltonian`: low-level hamiltonian object of type [`Hamiltonian`](@ref).
- `options`: options of the evolution in type [`KrylovOptions`](@ref).
"""
function CFETEvolution(reg, start_clock, durations, hamiltonian, options, algo)
    return CFETEvolution{typeof(reg),typeof(start_clock),typeof(hamiltonian),typeof(algo)}(
        reg,
        start_clock,
        durations,
        hamiltonian,
        options,
        algo,
    )
end

function Adapt.adapt_structure(to, x::CFETEvolution)
    return CFETEvolution(adapt(to, x.reg), x.start_clock, x.durations, adapt(to, x.hamiltonian), x.options, x.alg_table)
end

function CFETEvolution(reg::AbstractRegister, clocks, h, algo::CFETTables = CFET4_2(); kw...)
    all(≥(0), clocks) || throw(ArgumentError("clocks must not be negative"))
    options = from_kwargs(KrylovOptions; kw...)
    P = real(eltype(statevec(reg)))
    T = isreal(h) ? P : Complex{P}
    start_clock, durations = first(clocks), diff(clocks)
    
    ## checking if it is equal time:
    if length(unique(round.(durations,digits = 14) ) ) != 1
        throw(ArgumentError("durations must be equal (time slice must be equal)"))
    end 

    return CFETEvolution(reg, start_clock, durations, Hamiltonian(T, h, space(reg)), options, algo)
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


function emulate_step!(prob::CFETEvolution, step::Int, clock::Real, duration::Real)

    state = statevec(prob.reg)
    Ham = prob.hamiltonian
    
    ## each exponential-time prop. 
    for i in 1:length(prob.alg_table.Gs)
        
        #construct Ωi:
        Ωi = __construct_Ω(Ham, clock, duration, prob.alg_table, i)
        
        # perform evolution:
        prob.options.expmv_backend(-duration, im*Ωi, state; prob.options.tol)

    end

    # do we need this normalization? 
    if mod(step, prob.options.normalize_step) == 0
        normalize!(prob.reg)
    end

    if prob.options.normalize_finally && step == length(prob.durations)
        normalize!(prob.reg)
    end
    
    return prob

end

##==========================================================
#  Following are helper functions for calculating cfet tables
##==========================================================

function Lagendre(n::Int, x::Real)
    if n == 0
        return 1
    elseif n == 1
        return 2*x-1
    else
        return ((2n-1)*(2x-1)*Lagendre(n-1,x) - (n-1)*Lagendre(n-2,x))/n
    end
end

## arXiv:1102.5071 eq.(60)
#  Given Fs table, and Gaussian quadrature points xs and weights weights,
#  Generate Gs table with size (# of exponential time steps) x (# of quadrature points)
#  1. Fs have to be size (# of exponential time steps) x (# of order N/2+1)
#  2. (generally second axes of Fs it can be up to N but (R2) providing high order term can be ignore)

## TODO: calculate Lagendre on the fly
function __get_Gs(Fs::Vector{Vector{Float64}}, xs::Vector{Float64}, ws::Vector{Float64})
    
    Gs = Vector{Vector{Float64}}([zeros(length(xs)) for i in 1:length(Fs)])
    Nm = length(Fs[1])
    FL = length(Fs)

    ## Gs is filled in reverse order opposed from the definition becase eval of expm is in reverse order (right to left)
    for i in 1:FL
        for m in 1:length(xs)
            Gs[FL+1-i][m] = ws[m]*sum(Fs[i][n]*(2n-1)*Lagendre(n-1,xs[m]) for n in 1:Nm)
        end
    end
    return Gs

end
