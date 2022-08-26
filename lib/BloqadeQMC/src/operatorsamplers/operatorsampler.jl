abstract type AbstractOperatorSampler{K, T, P <: AbstractProbabilityVector{T}} end
firstindex(::AbstractOperatorSampler) = 1
lastindex(os::AbstractOperatorSampler) = length(os)
@inline normalization(os::AbstractOperatorSampler) = normalization(os.pvec)

struct OperatorSampler{K, T, P} <: AbstractOperatorSampler{K, T, P}
    operators::Vector{NTuple{K, Int}}
    pvec::P
    op_log_weights::Vector{T}
end


function OperatorSampler(operators::Vector{NTuple{K, Int}}, p::Vector{T}) where {K, T <: AbstractFloat}
    @assert length(operators) == length(p) "Given vectors must have the same length!"
    pvec = probability_vector(p)

    op_log_weights = log.(p)
    return OperatorSampler{K, T, typeof(pvec)}(operators, pvec, op_log_weights)
end

@inline rand(rng::AbstractRNG, os::OperatorSampler) = @inbounds os.operators[rand(rng, os.pvec)]

Base.@propagate_inbounds getlogweight(os::OperatorSampler, w::Int) = os.op_log_weights[w]

@inline length(os::OperatorSampler) = length(os.operators)

##############################################################################

