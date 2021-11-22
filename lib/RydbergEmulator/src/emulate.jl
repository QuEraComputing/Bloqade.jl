function emulate_routine!(r::RydbergReg{ComplexLayout}, t::Number, h::AbstractTerm, cache::AbstractMatrix)
    st = vec(r.state)
    update_term!(cache, h, r.subspace)
    expmv!(-im * t, cache, st)
    return r
end

function emulate_routine!(r::RydbergReg{RealLayout}, t::Number, h::AbstractTerm, cache::AbstractMatrix)
    error("discrete emulator does not support RealLayout yet")
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
struct DiscreteEmulationCache{Tv, Ti, S <: AbstractSparseMatrix{Tv, Ti}}
    H::S
end


# TODO: calculate the nnz colptr and rowval directly
"""
    DiscreteEmulationCache{Tv, Ti}(hs[, s::AbstractSpace=fullspace])

Create a `DiscreteEmulationCache`.

# Arguments

- `T`: element type of the storage.
- `hs`: a Hamiltonian expression term or a list of Hamiltonians.
- `s`: space type, default is [`fullspace`](@ref).
"""
function DiscreteEmulationCache{Tv, Ti}(h::AbstractTerm, s::AbstractSpace) where {Tv, Ti}
    is_time_dependent(h) && throw(ArgumentError("expect a time independent hamiltonian"))
    return DiscreteEmulationCache(SparseMatrixCSC{Tv, Ti}(h, s))
end

function DiscreteEmulationCache{Tv, Ti}(ts::AbstractVector, hs::Vector{<:AbstractTerm}, s::AbstractSpace) where {Tv, Ti}
    # use the one that has less zero term values
    _, idx = findmin(map(num_zero_term, hs))
    DiscreteEmulationCache{Tv, Ti}(hs[idx], s)
end

num_zero_term(t::Hamiltonian) = count(iszero, t.terms)
num_zero_term(t::AbstractTerm) = iszero(t) ? 1 : 0

"""
    emulate!(r, ts, hs)

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
"""
function emulate! end

struct PrecisionAdaptor{P} end

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

struct DiscreteEvolution{P, S, T <: Real, H <: AbstractTerm, Cache <: DiscreteEmulationCache}
    reg::S
    ts::Vector{T}
    hs::Vector{H}
    cache::Cache
    options::DiscreteOptions
end

Adapt.@adapt_structure DiscreteEvolution

get_space(r::Yao.ArrayReg) = fullspace
get_space(r::RydbergReg) = r.subspace

"""
    DiscreteEvolution{P}(register, ts, hs; kw...)

Create a `DiscreteEvolution` object that emulates a list of hamiltonians at discrete time steps
using Krylov subspace method, or trotterize a continuous function with `dt` then run the
trotterize integrator on it.

# Arguments

- `P`: optional, a type parameter that sets the problem precision type, default is
    the same as the `Yao.datatype` of given `register`.
- `register`: required, the evolution problem register, can be a [`RydbergReg`](@ref) or an `ArrayReg`
    from `Yao`.
- `ts`: required, the evolution time, can be a real number of a list of real numbers.
- `hs`: required, the evolution hamiltonian, can be one hamiltonian with time dependent parameters
    to trotterize or with constant parameters to evolve, or a list of hamiltonians with constant parameters.

# Keyword Arguments

- `cache`: discrete solver cache, see also [`DiscreteEmulationCache`](@ref).
- `normalize_step::Int`: run normalization per `normalize_step`, default is `5`.
- `normalize_finally::Bool`: normalize the state after the entire emulation ends, default is `true`.
- `progress::Bool`: show progress bar, default is `false`.
- `progress_step::Int`: update the progress bar per `progress_step`, default is `1`.
- `progress_name::String`: the printed name on progress bar, default is `"emulating"`.
- `dt::Real`: the time step of trotterization if `ts` is specified as
    a `Real` number and `hs` is a time dependent hamiltonian.
"""
function DiscreteEvolution{P}(
    r::Yao.AbstractRegister, ts::Vector{<:Real}, hs::Vector{<:AbstractTerm};
    index_type::Type=Cint,
    cache::DiscreteEmulationCache=default_discrete_evolution_cache(P, index_type, r, ts, hs),
    kw...) where P

    length(ts) == length(hs) || throw(ArgumentError("length of time does not match hamiltonian"))

    # convert to user specified precision type
    state = adapt(PrecisionAdaptor(P), r)
    options = DiscreteOptions(;kw...)
    # convert units
    ts = default_unit(μs, ts)
    ts = map(P, ts)

    S, T, H, C = typeof(state), eltype(ts), eltype(hs), typeof(cache)
    return DiscreteEvolution{P, S, T, H, C}(state, ts, hs, cache, options)
end

function default_discrete_evolution_cache(::Type{P}, ::Type{index_type}, r, ts, hs) where {P, index_type}
    all_real = if hs isa AbstractTerm
        isreal(hs)
    else
        all(isreal, hs)
    end
    T = all_real ? P : Complex{P}
    return DiscreteEmulationCache{T, index_type}(ts, hs, get_space(r))
end

function DiscreteEvolution(r::Yao.AbstractRegister, ts, hs; kw...)
    return DiscreteEvolution{real(Yao.datatype(r))}(r, ts, hs; kw...)
end

# TODO: use GarishPrint after it gets smarter
function Base.show(io::IO, mime::MIME"text/plain", prob::DiscreteEvolution{P}) where P
    indent = get(io, :indent, 0)
    tab(indent) = " "^indent
    println(io, tab(indent), "DiscreteEvolution{", P, "}:")
    
    # state info
    print(io, tab(indent), "  reg: ")
    printstyled(io, typeof(prob.reg); color=:green)
    println(io)
    print(io, tab(indent), "  reg storage: ")
    printstyled(io, Base.format_bytes(sizeof(Yao.state(prob.reg))); color=:yellow)
    println(io)
    println(io)

    print(io, tab(indent), "  time: ")
    println(io, prob.ts, " μs")

    println(io)
    println(io, tab(indent), "  hamiltonian: ")
    if length(prob.hs) < 6
        for h in prob.hs
            show(IOContext(io, :indent=>indent+4), mime, h)
            println(io)
        end
    else
        show(IOContext(io, :indent=>indent+4), mime, first(prob.hs))
        println(io)
        println(io)
        println(io, tab(indent), "       ⋮")
        println(io)
        show(IOContext(io, :indent=>indent+4), mime, last(prob.hs))
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

storage_size(r::RydbergReg) = sizeof(r.state) + sizeof(r.subspace.subspace_v) + sizeof(r.subspace.map)
storage_size(r::Yao.ArrayReg) = sizeof(r.state)
storage_size(H::Array) = sizeof(H)
storage_size(x) = sizeof(x) # fallback to sizeof

Base.@propagate_inbounds function emulate_step!(
    prob::DiscreteEvolution,
    ts::Vector{<:Number}, hs::Vector{<:AbstractTerm},
    i::Int)

    emulate_routine!(prob.reg, ts[i], hs[i], prob.cache.H)

    if mod(i, prob.options.normalize_step) == 0
        normalize!(prob.reg)
    end
    return prob
end

# constant params
function emulate_step!(
    prob::DiscreteEvolution, t::Number, h::AbstractTerm, ::Int)
    # i is a constant 1
    emulate_routine!(prob.reg, t, h, prob.cache.H)    
    return prob
end

# continuous function
Base.@propagate_inbounds function emulate_step!(
    prob::DiscreteEvolution,
    t::AbstractVector{<:Number}, h::AbstractTerm,
    i::Int)

    emulate_routine!(prob.reg, step(t), h(t[i]), prob.cache.H)
    if mod(i, prob.options.normalize_step) == 0
        normalize!(prob.reg)
    end
    return prob
end

function emulate!(prob::DiscreteEvolution)
    niterations = length(prob.ts)
    @inbounds if prob.options.progress
        ProgressLogging.progress() do id
            for idx in 1:niterations
                emulate_step!(prob, prob.ts, prob.hs, idx)
                if prob.options.progress && mod(idx, prob.options.progress_step) == 0
                    @info prob.options.progress_name progress=idx/niterations _id=id
                end
            end
        end
    else
        for idx in 1:niterations
            emulate_step!(prob, prob.ts, prob.hs, idx)
        end
    end

    if prob.options.normalize_finally
        normalize!(prob.reg)
    end
    return prob
end

function emulate!(r::Yao.AbstractRegister, ts::Vector{<:Number}, hs::Vector{<:AbstractTerm}; kw...)
    prob = DiscreteEvolution(r, ts, hs; kw...)
    emulate!(prob)
    return prob.reg
end

"""
    emulate(ts, hs; kwargs...)

Non in-place version of [`emulate!`](@ref). See [`emulate!`](@ref) for valid kwargs.
"""
function emulate(ts, hs; kwargs...)
    precision_t = real(eltype(ts))
    r = Yao.zero_state(Complex{precision_t}, nsites(hs))
    return emulate!(r, ts, hs; kwargs...)
end

"""
    emulate(s::Subspace, ts, hs; kwargs...)

Non in-place version of [`emulate!`](@ref). See [`emulate!`](@ref) for valid kwargs.
"""
function emulate(s::Subspace, ts, hs; kwargs...)
    @assert s.nqubits == nsites(hs) "qubit mismatch"
    precision_t = real(eltype(ts))
    r = zero_state(Complex{precision_t}, s)
    return emulate!(r, ts, hs; kwargs...)
end
