function emulate_routine!(r::RydbergReg, t::Number, h::AbstractTerm, cache::AbstractMatrix)
    st = vec(r.state)
    update_term!(cache, h, r.subspace)
    expmv!(-im * t, cache, st)
    return r
end

function emulate_routine!(r::Yao.ArrayReg, t::Number, h::AbstractTerm, cache::AbstractMatrix)
    st = vec(r.state)
    update_term!(cache, h)
    expmv!(-im * t, cache, st)
    return r
end

"""
    DiscreteEmulationCache{C}

Cache type for the discrete emulation, type variable `C` is the actual
matrix type for the hamiltonian storage.
    
When we do Krylov subspace based time evolution emulation with a sequence
of hamiltonians of similar terms, e.g all the hamiltonians are
`RydInteract + XTerm + ZTerm`, the sparse structure will be the same, thus
we can re-use the sparse matrix generated for the first hamiltonian again via
[`update_term!`](@ref) in the following calculation to reduce the memory usage
and speed up the emulation.
"""
struct DiscreteEmulationCache{C}
    H::C
end


# TODO: calculate the nnz colptr and rowval directly
"""
    DiscreteEmulationCache(::Type{T}, h_or_hs[, s::AbstractSpace=fullspace])

Create a `DiscreteEmulationCache`.

# Arguments

- `T`: element type of the storage.
- `h_or_hs`: a Hamiltonian expression term or a list of Hamiltonians.
- `s`: space type, default is [`fullspace`](@ref).
"""
function DiscreteEmulationCache(::Type{T}, h::AbstractTerm, s::AbstractSpace=fullspace) where {T}
    is_time_dependent(h) && throw(ArgumentError("expect a time independent hamiltonian"))
    return DiscreteEmulationCache(SparseMatrixCSC{T}(h, s))
end

function DiscreteEmulationCache(::Type{T}, ts::AbstractVector, hs::Vector{<:AbstractTerm}, s::AbstractSpace=fullspace) where {T}
    # use the one that has less zero term values
    _, idx = findmin(map(num_zero_term, hs))
    DiscreteEmulationCache(T, hs[idx], s)
end

function DiscreteEmulationCache(::Type{T}, ts::AbstractVector, h::AbstractTerm, s::AbstractSpace=fullspace) where {T}
    # use the one that has less zero term values
    _, idx = findmin(map(t->num_zero_term(h(t)), ts))
    DiscreteEmulationCache(T, h(ts[idx]), s)
end

function DiscreteEmulationCache(::Type{T}, t::Number, h::AbstractTerm, s::AbstractSpace=fullspace) where {T}
    DiscreteEmulationCache(T, h, s)
end

num_zero_term(t::Hamiltonian) = count(iszero, t.terms)
num_zero_term(t::AbstractTerm) = iszero(t) ? 1 : 0

"""
    emulate!(r, ts, hs[; cache=DiscreteEmulationCache(ts, hs)])

Emulate the time evolution of a sequence of Hamiltonians `hs` of time length `ts` for each Hamiltonian.
An optional argument can be feeded to preallocate the memory.

    emulate!(r, t::Real, h::AbstractTerm[; algo=Tsit5(), cache=<default cache>, kwargs...])

Emulate the continuous time evolution of given Hamiltonian `h` that has dependency on `t`. The default algorithm
we use is `Tsit5`. Available `kwargs` can be found in [Common Solver Options](https://diffeq.sciml.ai/latest/basics/common_solver_opts/).
Available algorithms can be found at [Full List of Methods](https://diffeq.sciml.ai/latest/solvers/ode_solve/#Full-List-of-Methods).

# Emulation Cache

Emulation of time evolution requires allocating intermediate variables that can be large. Thus preallocating
these intermediate variables and share this memory between iterations can speed up the simulation when `emulate!`
is called for multiple times for similar Hamiltonian. By similar Hamiltonian, we mean Hamiltonians that contains
the same term but can have different parameters.

## Time Independent Cache

The cache for time independent simulation contains two part: a) a `SparseMatrixCSC` matrix to store the intermediate
Hamiltonian matrices. b) a `KrylovSubspace` object to store the Krylov subspace.

## continuous Time Cache

For continuous time emulation, one only needs to create a `SparseMatrixCSC` matrix to store the intermediate
Hamiltonian matrices.

## Examples

The simplest way to use cache inside a large loop is via closure. The following example
returns a closure that takes time `ts` as parameters during simulation.

```jl
function your_emulation_task(n, subspace, hs)
    r = zero_state(n, subspace)
    cache = DiscreteEmulationCache(first(hs), subspace)
    return function task(ts)
        set_zero_state!(r)
        emulate!(r, ts, hs; cache=cache)
    end
end
```

or if you just want to use `ts` and `ϕs` as your parameters

```jl
function your_emulation_task(n, subspace)
    r = zero_state(n, subspace)
    cache = DiscreteEmulationCache(simple_rydberg(n, 1.0), subspace)

    return function task(xs)
        set_zero_state!(r)
        ts = xs[1:length(xs)÷2]
        ϕs = xs[length(xs)÷2+1:end]
        hs = simple_rydberg.(n, ϕs)
        emulate!(r, ts, hs; cache=cache)
    end
end
```
"""
function emulate! end

struct PrecisionAdaptor{P <: AbstractFloat} end

PrecisionAdaptor(P::Type) = PrecisionAdaptor{P}()
Adapt.adapt_storage(::PrecisionAdaptor{P}, x::Real) where P = P(x)
Adapt.adapt_storage(::PrecisionAdaptor{P}, x::Complex) where P = Complex{P}(x)
Adapt.adapt_storage(::PrecisionAdaptor{P}, x::Array) where P = convert(Array{P}, x)
Adapt.adapt_storage(::PrecisionAdaptor{P}, x::Array{<:Complex}) where P = convert(Array{Complex{P}}, x)

abstract type EmulationOptions end

@option struct DiscreteOptions <: EmulationOptions
    progress::Bool = false
    progress_step::Int = 1
    progress_name::String = "emulating"
    normalize_step::Int = 5
    normalize_finally::Bool = true
end

struct DiscreteEvolution{P, S, T, H, C}
    state::S
    t_or_ts::T
    h_or_hs::H
    cache::DiscreteEmulationCache{C}
    options::DiscreteOptions
end

Adapt.@adapt_structure DiscreteEvolution

get_space(r::Yao.ArrayReg) = fullspace
get_space(r::RydbergReg) = r.subspace

"""
    DiscreteEvolution{P}(register, t_or_ts, h_or_hs; kw...)

Create a `DiscreteEvolution` object that emulates a list of hamiltonians at discrete time steps
using Krylov subspace method, or trotterize a continuous function with `dt` then run the
trotterize integrator on it.

# Arguments

- `P`: optional, a type parameter that sets the problem precision type, default is
    the same as the `Yao.datatype` of given `register`.
- `register`: required, the evolution problem register, can be a [`RydbergReg`](@ref) or an `ArrayReg`
    from `Yao`.
- `t_or_ts`: required, the evolution time, can be a real number of a list of real numbers.
- `h_or_hs`: required, the evolution hamiltonian, can be one hamiltonian with time dependent parameters
    to trotterize or with constant parameters to evolve, or a list of hamiltonians with constant parameters.

# Keyword Arguments

- `cache`: discrete solver cache, see also [`DiscreteEmulationCache`](@ref).
- `normalize_step::Int`: run normalization per `normalize_step`, default is `5`.
- `normalize_finally::Bool`: normalize the state after the entire emulation ends, default is `true`.
- `progress::Bool`: show progress bar, default is `true`.
- `progress_step::Int`: update the progress bar per `progress_step`, default is `1`.
- `progress_name::String`: the printed name on progress bar, default is `"emulating"`.
- `dt::Real`: the time step of trotterization if `t_or_ts` is specified as
    a `Real` number and `h_or_hs` is a time dependent hamiltonian.
"""
function DiscreteEvolution{P}(
    r::Yao.AbstractRegister, t_or_ts, h_or_hs;
    cache=nothing, dt::Real=1e-3, kw...) where P

    if !(h_or_hs isa AbstractTerm)
        length(t_or_ts) == length(h_or_hs) || throw(ArgumentError("length of time does not match hamiltonian"))
    end

    eltype(t_or_ts) <: Real || throw(ArgumentError("time should be a real number, got $(eltype(t_or_ts))"))

    state = adapt(PrecisionAdaptor(P), r)
    options = DiscreteOptions(;kw...)

    # always convert it to a range if h_or_hs is time dependent
    if t_or_ts isa Real && is_time_dependent(h_or_hs)
        t_or_ts = zero(t_or_ts):dt:t_or_ts
    end

    # convert units
    t_or_ts = default_unit(μs, t_or_ts)
    # convert to given precision
    t_or_ts = map(P, t_or_ts)

    if isnothing(cache)
        all_real = if h_or_hs isa AbstractTerm
            isreal(h_or_hs)
        else
            all(isreal, h_or_hs)
        end
        T = all_real ? P : Complex{P}
        cache = DiscreteEmulationCache(T, t_or_ts, h_or_hs, get_space(r))
    end

    S, T, H, C = typeof(state), typeof(t_or_ts), typeof(h_or_hs), typeof(cache.H)
    return DiscreteEvolution{P, S, T, H, C}(state, t_or_ts, h_or_hs, cache, options)
end

function DiscreteEvolution(r::Yao.AbstractRegister, t_or_ts, h_or_hs; kw...)
    return DiscreteEvolution{real(Yao.datatype(r))}(r, t_or_ts, h_or_hs; kw...)
end

# TODO: use GarishPrint after it gets smarter
function Base.show(io::IO, mime::MIME"text/plain", prob::DiscreteEvolution{P}) where P
    indent = get(io, :indent, 0)
    tab(indent) = " "^indent
    println(io, tab(indent), "DiscreteEvolution{", P, "}:")
    
    # state info
    print(io, tab(indent), "  state: ")
    printstyled(io, typeof(prob.state); color=:green)
    println(io)
    print(io, tab(indent), "  state storage: ")
    printstyled(io, Base.format_bytes(sizeof(Yao.state(prob.state))); color=:yellow)
    println(io)
    println(io)

    print(io, tab(indent), "  time: ")
    if prob.t_or_ts isa AbstractRange
        printstyled(io, typeof(prob.t_or_ts); color=:green)
        println(io)
        println(io, tab(indent), "    start: ", first(prob.t_or_ts), " μs")
        println(io, tab(indent), "     step: ", step(prob.t_or_ts), " μs")
        println(io, tab(indent), "     stop: ", last(prob.t_or_ts), " μs")
    else
        println(io, prob.t_or_ts, " μs")
    end

    println(io)
    println(io, tab(indent), "  hamiltonian: ")
    if prob.h_or_hs isa AbstractTerm
        show(IOContext(io, :indent=>indent+4), mime, prob.h_or_hs)
    elseif length(prob.h_or_hs) < 6
        for h in prob.h_or_hs
            show(IOContext(io, :indent=>indent+4), mime, h)
            println(io)
        end
    else
        show(IOContext(io, :indent=>indent+4), mime, first(prob.h_or_hs))
        println(io)
        println(io)
        println(io, tab(indent), "       ⋮")
        println(io)
        show(IOContext(io, :indent=>indent+4), mime, last(prob.h_or_hs))
    end
    println(io)
    println(io)
    print(io, tab(indent), "  hamiltonian storage: ")
    printstyled(io, Base.format_bytes(storage_size(prob.cache)); color=:yellow)
    println(io)

    println(io)
    println(io, tab(indent), "  options: ")

    nfs = nfields(prob.options)
    for idx in 1:nfs
        name = fieldname(DiscreteOptions, idx)
        print(io, tab(indent), "    $name: ", repr(getfield(prob.options, idx)))
        if idx != nfs
            println(io)
        end
    end
end

function storage_size(cache::DiscreteEmulationCache)
    return storage_size(cache.H)
end

function storage_size(H::SparseMatrixCSC)
    return sizeof(H.colptr) + sizeof(H.rowval) + sizeof(H.nzval)
end

storage_size(H::Array) = sizeof(H)

Base.@propagate_inbounds function emulate_step!(
    prob::DiscreteEvolution,
    t_or_ts::Vector{<:Number}, h_or_hs::Vector{<:AbstractTerm},
    i::Int)

    emulate_routine!(prob.state, t_or_ts[i], h_or_hs[i], prob.cache.H)

    if mod(i, prob.options.normalize_step) == 0
        normalize!(prob.state)
    end
    return prob
end

# constant params
function emulate_step!(
    prob::DiscreteEvolution, t::Number, h::AbstractTerm, ::Int)
    # i is a constant 1
    emulate_routine!(prob.state, t, h, prob.cache.H)    
    return prob
end

# continuous function
Base.@propagate_inbounds function emulate_step!(
    prob::DiscreteEvolution,
    t::AbstractVector{<:Number}, h::AbstractTerm,
    i::Int)

    emulate_routine!(prob.state, step(t), h(t[i]), prob.cache.H)
    if mod(i, prob.options.normalize_step) == 0
        normalize!(prob.state)
    end
    return prob
end

function emulate!(prob::DiscreteEvolution)
    niterations = length(prob.t_or_ts)
    @inbounds if prob.options.progress
        ProgressLogging.progress() do id
            for idx in 1:niterations
                emulate_step!(prob, prob.t_or_ts, prob.h_or_hs, idx)
                if prob.options.progress && mod(idx, prob.options.progress_step) == 0
                    @info prob.options.progress_name progress=idx/niterations _id=id
                end
            end
        end
    else
        for idx in 1:niterations
            emulate_step!(prob, prob.t_or_ts, prob.h_or_hs, idx)
        end
    end

    if prob.options.normalize_finally
        normalize!(prob.state)
    end
    return prob
end

function emulate!(r::Yao.AbstractRegister, ts::Vector{<:Number}, hs::Vector{<:AbstractTerm}; kw...)
    prob = DiscreteEvolution(r, ts, hs; kw...)
    emulate!(prob)
    return prob.state
end

"""
    emulate(t_or_ts, h_or_hs; kwargs...)

Non in-place version of [`emulate!`](@ref). See [`emulate!`](@ref) for valid kwargs.
"""
function emulate(t_or_ts, h_or_hs; kwargs...)
    precision_t = real(eltype(t_or_ts))
    r = Yao.zero_state(Complex{precision_t}, nsites(h_or_hs))
    return emulate!(r, t_or_ts, h_or_hs; kwargs...)
end

"""
    emulate(s::Subspace, t_or_ts, h_or_hs; kwargs...)

Non in-place version of [`emulate!`](@ref). See [`emulate!`](@ref) for valid kwargs.
"""
function emulate(s::Subspace, t_or_ts, h_or_hs; kwargs...)
    precision_t = real(eltype(t_or_ts))
    r = zero_state(Complex{precision_t}, nsites(h_or_hs), s)
    return emulate!(r, t_or_ts, h_or_hs; kwargs...)
end
