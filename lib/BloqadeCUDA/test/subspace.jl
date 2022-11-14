using Test
using Random
using CUDA
using YaoSubspaceArrayReg
using YaoArrayRegister
using StatsBase
using CuYao
using GPUArrays
using Adapt
using BloqadeCUDA


struct Sampler{A,B}
    cum_prob::A
    values::B
end

struct CuWeights{T<:Real}
    values::CuVector{T}
    sum::T
end


function Adapt.adapt_structure(to, sampler::Sampler)
    cum_prob = Adapt.adapt_structure(to, sampler.cum_prob)
    values = Adapt.adapt_structure(to, sampler.values)
    Sampler(cum_prob,values)
end


function StatsBase.Weights(values::CuVector{T}) where {T<:Real}
    return CuWeights{eltype(values)}(values,sum(values))
end

function (sampler::Sampler)(x)
    i = searchsortedfirst(sampler.cum_prob, x)
    i = clamp(i, firstindex(sampler.values), lastindex(sampler.values))
    @inbounds sampler.values[i]
end


function sample(rng::AbstractRNG, subspace_v::CuVector, weights::CuWeights,nshots::Integer)
    dices = rand(rng, nshots)
    sampler = Sampler(cumsum(weights.values),subspace_v)

    Array(sampler.(dices))
end



space = Subspace(10, sort(randperm(1 << 10)[1:76] .- 1))
r = SubspaceArrayReg(randn(ComplexF64, 76), space)
normalize!(r)
dr = cu(r)


weights = Weights(abs2.(relaxedvec(dr)))
subspace_v = vec(dr.subspace)

samples = sample(CURAND.default_rng(),subspace_v,weights,1000)

# weights = Weights(abs2.(relaxedvec(dr)))
# @which measure(r)

# measure(dr)