abstract type AbstractProbabilityVector{T <: Real} <: AbstractVector{T} end
Base.@propagate_inbounds getindex(pvec::AbstractProbabilityVector, i) = getweight(pvec, i)
Base.@propagate_inbounds setindex!(pvec::AbstractProbabilityVector, w, i) = setweight!(pvec, w, i)
Base.@propagate_inbounds setprobability!(pvec::AbstractProbabilityVector{T}, p::T, i::Int) where T = setweight!(pvec, p*normalization(pvec), i)
Base.@propagate_inbounds getprobability(pvec::AbstractProbabilityVector, i) = getweight(pvec, i) / normalization(pvec)

firstindex(::AbstractProbabilityVector) = 1
lastindex(pvec::AbstractProbabilityVector) = length(pvec)
size(pvec::AbstractProbabilityVector) = (length(pvec),)

###############################################################################

include("alias.jl")

###############################################################################

# Get the "recommended" probability vector type
probability_vector(p::Vector{T}) where T = ProbabilityAlias(p)
