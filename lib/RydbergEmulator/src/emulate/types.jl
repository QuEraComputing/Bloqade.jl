"""
    KrylovEmulationCache{C}

Cache type for the discrete emulation, type variable `C` is the actual
matrix type for the hamiltonian storage.
    
When we do Krylov subspace based time evolution emulation with a sequence
of hamiltonians of similar terms, e.g all the hamiltonians are
`RydInteract + XTerm + ZTerm`, the sparse structure will be the same, thus
we can re-use the sparse matrix generated for the first hamiltonian again via
[`update_term!`](@ref) in the following calculation to reduce the memory usage
and speed up the emulation.
"""
struct KrylovEmulationCache{Tv, Ti, S <: AbstractSparseMatrix{Tv, Ti}}
    H::S
end


# TODO: calculate the nnz colptr and rowval directly
"""
    KrylovEmulationCache{Tv, Ti}(hs[, s::AbstractSpace=fullspace])

Create a `KrylovEmulationCache`.

# Arguments

- `T`: element type of the storage.
- `hs`: a Hamiltonian expression term or a list of Hamiltonians.
- `s`: space type, default is [`fullspace`](@ref).
"""
function KrylovEmulationCache{Tv, Ti}(h::AbstractTerm, s::AbstractSpace) where {Tv, Ti}
    is_time_dependent(h) && throw(ArgumentError("expect a time independent hamiltonian"))
    return KrylovEmulationCache(SparseMatrixCSC{Tv, Ti}(h, s))
end

function KrylovEmulationCache{Tv, Ti}(hs::Vector{<:AbstractTerm}, s::AbstractSpace) where {Tv, Ti}
    # use the one that has less zero term values
    _, idx = findmin(map(num_zero_term, hs))
    KrylovEmulationCache{Tv, Ti}(hs[idx], s)
end

abstract type EmulationOptions end

@option struct DiscreteOptions <: EmulationOptions
    progress::Bool = false
    progress_step::Int = 1
    progress_name::String = "emulating"
    normalize_step::Int = 5
    normalize_finally::Bool = true
end

struct KrylovEvolution{P, S, T <: Real, H <: AbstractTerm, Cache <: KrylovEmulationCache}
    reg::S
    durations::Vector{T}
    hs::Vector{H}
    cache::Cache
    options::DiscreteOptions
end

Adapt.@adapt_structure KrylovEvolution

"""
    KrylovEvolution{P}(register, durations, hs; kw...)

Create a `KrylovEvolution` object that emulates a list of hamiltonians at discrete time steps
using Krylov subspace method, or trotterize a continuous function with `dt` then run the
trotterize integrator on it.

# Arguments

- `P`: optional, a type parameter that sets the problem precision type, default is
    the same as the `Yao.datatype` of given `register`.
- `register`: required, the evolution problem register, can be a [`RydbergReg`](@ref) or an `ArrayReg`
    from `Yao`.
- `durations`: required, the evolution durations of each hamiltonian, should be a list of real numbers.
- `hs`: required, the evolution hamiltonian, a list of hamiltonians with constant parameters.

# Keyword Arguments

- `cache`: discrete solver cache, see also [`KrylovEmulationCache`](@ref).
- `normalize_step::Int`: run normalization per `normalize_step`, default is `5`.
- `normalize_finally::Bool`: normalize the state after the entire emulation ends, default is `true`.
- `progress::Bool`: show progress bar, default is `false`.
- `progress_step::Int`: update the progress bar per `progress_step`, default is `1`.
- `progress_name::String`: the printed name on progress bar, default is `"emulating"`.
- `dt::Real`: the time step of trotterization if `ts` is specified as
    a `Real` number and `hs` is a time dependent hamiltonian.
"""
function KrylovEvolution{P}(
    r::Yao.AbstractRegister, durations::Vector{<:Real}, hs::Vector{<:AbstractTerm};
    index_type::Type=Cint,
    cache::KrylovEmulationCache=default_discrete_evolution_cache(P, index_type, r, durations, hs),
    kw...) where P

    all(hs) do h
        nsites(h) == nsites(first(hs))
    end || throw(ArgumentError("number of sites does not match"))
    length(durations) == length(hs) || throw(ArgumentError("length of time does not match hamiltonian"))

    # convert to user specified precision type
    state = adapt(PrecisionAdaptor(P), r)
    options = DiscreteOptions(;kw...)
    # convert units
    durations = default_unit(μs, durations)
    durations = map(P, durations)

    S, T, H, C = typeof(state), eltype(durations), eltype(hs), typeof(cache)
    return KrylovEvolution{P, S, T, H, C}(state, durations, hs, cache, options)
end

function KrylovEvolution(durations::Vector{<:Real}, hs::Vector{<:AbstractTerm}; kw...)
    return KrylovEvolution(zero_state(Complex{eltype(durations)}, nsites(first(hs))), durations, hs; kw...)
end

function KrylovEvolution(r::Yao.AbstractRegister, durations::Vector{<:Real}, hs::Vector{<:AbstractTerm}; kw...)
    return KrylovEvolution{real(Yao.datatype(r))}(r, durations, hs; kw...)
end

function default_discrete_evolution_cache(::Type{P}, ::Type{index_type}, r, durations, hs) where {P, index_type}
    all_real = if hs isa AbstractTerm
        isreal(hs)
    else
        all(isreal, hs)
    end
    T = all_real ? P : Complex{P}
    return KrylovEmulationCache{T, index_type}(hs, get_space(r))
end
