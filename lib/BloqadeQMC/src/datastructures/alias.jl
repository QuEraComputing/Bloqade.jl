###############################################################################
# Walker's Alias Method: draw samples in O(1) time, initialize vector in O(n)
#   caveat: unsure if it's possible to build an efficient update scheme
#           for these types of probability vectors

probability_vector(p::Vector{T}) where T = ProbabilityAlias(p)

struct ProbabilityAlias{T} # <: AbstractProbabilityVector{T}     We should be fine without abstract subtyping   
    normalization::T

    cutoffs::Vector{T}
    alias::Vector{Int}

    # initialize using Vose's algorithm
    function ProbabilityAlias{T}(weights::Vector{T}) where {T <: AbstractFloat}
        if length(weights) == 0
            throw(ArgumentError("weight vector must have non-zero length!"))
        end
        if any(x -> x < zero(T), weights)
            throw(ArgumentError("weights must be non-negative!"))
        end

        weights_ = copy(weights)
        normalization = sum(weights)
        N = length(weights)
        avg = normalization / N

        underfull = Stack{Int}()
        overfull = Stack{Int}()

        for (i, w) in enumerate(weights_)
            if w >= avg
                push!(overfull, i)
            else
                push!(underfull, i)
            end
        end

        cutoffs = zeros(float(T), N)
        alias = zeros(Int, N)

        while !isempty(underfull) && !isempty(overfull)
            less = pop!(underfull)
            more = pop!(overfull)

            cutoffs[less] = weights_[less] * N
            alias[less] = more

            weights_[more] += weights_[less]
            weights_[more] -= avg

            if weights_[more] >= avg
                push!(overfull, more)
            else
                push!(underfull, more)
            end
        end

        while !isempty(underfull)
            cutoffs[pop!(underfull)] = normalization
        end

        while !isempty(overfull)
            cutoffs[pop!(overfull)] = normalization
        end

        cutoffs /= normalization

        new{T}(sum(weights), cutoffs, alias)
    end
end

ProbabilityAlias(p::Vector{T}) where T = ProbabilityAlias{T}(p)
@inline length(pvec::ProbabilityAlias) = length(pvec.cutoffs)
@inline normalization(pvec::ProbabilityAlias) = pvec.normalization

# For safety, copying the following over from probabilityvector.jl. Probably not needed:

# firstindex(::AbstractProbabilityVector) = 1
# lastindex(pvec::AbstractProbabilityVector) = length(pvec)
# size(pvec::AbstractProbabilityVector) = (length(pvec),)

function show(io::IO, p::ProbabilityAlias{T}) where T
    r = repr(p.normalization; context=IOContext(io, :limit=>true))
    print(io, "ProbabilityAlias{$T}($r)")
end

# function Base.rand(rng::AbstractRNG, pvec::ProbabilityAlias{T}) where T
#     u, i::Int = modf(muladd(length(pvec), rand(rng), 1.0))
#     return @inbounds (u < pvec.cutoffs[i]) ? i : pvec.alias[i]
# end

# xorshift prngs seem to be noticably faster if you just sample from them twice
function Base.rand(rng::AbstractRNG, pvec::ProbabilityAlias{T}) where T
    u = rand(rng)
    i = rand(rng, 1:length(pvec))
    return @inbounds (u < pvec.cutoffs[i]) ? i : pvec.alias[i]
end
