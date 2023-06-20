"""
    struct KrylovEvolution
    KrylovEvolution(reg::AbstractRegister, clocks, h; kw...)

Create a `KrylovEvolution` object that describes a time evolution
using Krylov subspace methods.

# Arguments

- `reg`: a register, should be a subtype of `AbstractRegister`.
- `clocks`: the clocks of this time evolution at each step.
- `h`: a hamiltonian expression.

# Keyword Arguments

- `progress`: show progress bar, default is `false`.
- `progress_name`: progress bar name, default is `"emulating"`.
- `normalize_step`: normalize the state every `normalize_step`.
- `normalize_finally`: wether normalize the state in the end of evolution, default is `true`.
- `tol`: tolerance of the Krylov method, default is `1e-7`

# Examples

The following is the simplest way of using `KrylovEvolution`
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


julia> prob = KrylovEvolution(r, 0.0:1e-2:0.1, h);

julia> emulate!(prob); # run the emulation
```
"""
struct KrylovEvolution{Reg<:AbstractRegister,T<:Real,H<:Hamiltonian} <: Evolver
    reg::Reg
    start_clock::T
    durations::Vector{T}
    hamiltonian::H
    options::KrylovOptions

    function KrylovEvolution{Reg,T,H}(reg, start_clock, durations, hamiltonian, options) where {Reg,T,H}
        start_clock ≥ 0 || throw(ArgumentError("start clock must not be negative"))
        all(≥(0), durations) || throw(ArgumentError("durations must not be negative"))
        return new{Reg,T,H}(reg, start_clock, durations, hamiltonian, options)
    end
end

"""
    KrylovEvolution(reg, start_clock, durations, hamiltonian, options)

Create a `KrylovEvolution` object.

# Arguments

- `reg`: a register object.
- `start_clock`: start clock of the evolution.
- `durations`: list of durations at each time step.
- `hamiltonian`: low-level hamiltonian object of type [`Hamiltonian`](@ref).
- `options`: options of the evolution in type [`KrylovOptions`](@ref).
"""
function KrylovEvolution(reg, start_clock, durations, hamiltonian, options)
    return KrylovEvolution{typeof(reg),typeof(start_clock),typeof(hamiltonian)}(
        reg,
        start_clock,
        durations,
        hamiltonian,
        options,
    )
end

function Adapt.adapt_structure(to, x::KrylovEvolution)
    return KrylovEvolution(adapt(to, x.reg), x.start_clock, x.durations, adapt(to, x.hamiltonian), x.options)
end

function KrylovEvolution(reg::AbstractRegister, clocks, h; kw...)
    all(≥(0), clocks) || throw(ArgumentError("clocks must not be negative"))
    options = from_kwargs(KrylovOptions; kw...)
    P = real(eltype(statevec(reg)))
    T = isreal(h) ? P : Complex{P}
    start_clock, durations = first(clocks), diff(clocks)
    return KrylovEvolution(reg, start_clock, durations, Hamiltonian(T, h, space(reg)), options)
end

function emulate_step!(prob::KrylovEvolution, step::Int, clock::Real, duration::Real)
    state = statevec(prob.reg)
    h = prob.hamiltonian

    expmv!(-duration * im, h(clock), state; prob.options.tol)
    if mod(step, prob.options.normalize_step) == 0
        normalize!(prob.reg)
    end

    if prob.options.normalize_finally && step == length(prob.durations)
        normalize!(prob.reg)
    end
    return prob
end





