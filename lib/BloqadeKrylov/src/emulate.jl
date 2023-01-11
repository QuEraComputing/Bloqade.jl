"""
    struct KrylovOptions

Krylov evolution options.
"""
@option struct KrylovOptions
    progress::Bool = false
    progress_step::Int = 1
    progress_name::String = "emulating"
    normalize_step::Int = 5
    normalize_finally::Bool = true
    tol::Float64 = 1e-7
end

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
struct KrylovEvolution{Reg<:AbstractRegister,T<:Real,H<:Hamiltonian}
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

Base.length(prob::KrylovEvolution) = length(prob.durations) + 1

function Base.iterate(prob::KrylovEvolution)
    info = (; step = 1, reg = prob.reg, clock = prob.start_clock, duration = zero(prob.start_clock))
    return info, (2, prob.start_clock)
end

Base.@propagate_inbounds function Base.iterate(prob::KrylovEvolution, (step, clock))
    step > length(prob) && return

    duration = prob.durations[step-1]
    emulate_step!(prob, step, clock, duration)

    info = (; step, reg = prob.reg, clock = clock + duration, duration)
    return info, (step + 1, clock + duration)
end

function BloqadeExpr.emulate!(prob::KrylovEvolution)
    niterations = length(prob)
    @inbounds if prob.options.progress
        ProgressLogging.progress() do id
            for info in prob
                if prob.options.progress && mod(info.step, prob.options.progress_step) == 0
                    @info prob.options.progress_name progress = info.step / niterations _id = id
                end
            end
        end
    else
        for info in prob
        end
    end
    return prob
end

tab(indent) = " "^indent

function Base.show(io::IO, mime::MIME"text/plain", prob::KrylovEvolution)
    indent = get(io, :indent, 0)
    println(io, tab(indent), "KrylovEvolution:")
    # state info
    print_state_info(io, prob)
    println(io)

    # clocks
    println(io, tab(indent + 2), "clocks")
    println(io, tab(indent + 4), "start:", prob.start_clock, "μs")
    println(io, tab(indent + 4), " last:", prob.start_clock + sum(prob.durations), "μs")
    println(io)

    # equation info
    show(IOContext(io, :indent => indent + 2), mime, prob.hamiltonian)
    println(io)
    println(io)

    println(io, tab(indent + 2), "Options:")
    for name in fieldnames(KrylovOptions)
        println(io, tab(indent + 4), name, "=", repr(getfield(prob.options, name)))
    end
end

function print_state_info(io::IO, prob::KrylovEvolution)
    indent = get(io, :indent, 0)
    println(io, tab(indent + 2), "register info:")
    print(io, tab(indent + 4), "type: ")
    printstyled(io, typeof(prob.reg); color = :green)
    println(io)

    print(io, tab(indent + 4), "storage size: ")
    printstyled(io, Base.format_bytes(storage_size(prob.reg)); color = :yellow)
    return println(io)
end
