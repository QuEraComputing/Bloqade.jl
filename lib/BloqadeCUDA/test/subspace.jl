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

# TODO: use Alias sampler
struct SamplerWithValues{T,A<:AbstractVector{T},B<:AbstractVector}
    total::T
    cumulative_weights::A
    values::B
end

function SamplerWithValues(weights,values)
    length(weights) != length(values) && error("Sampler expects values with same length as weights.")
    total = sum(weights)
    cumulative_weights = cumsum(weights)
    return SamplerWithValues(total,cumulative_weights,values)
end

function Adapt.adapt_structure(to, sampler::SamplerWithValues)
    cumulative_weights = Adapt.adapt_structure(to, sampler.cumulative_weights)
    values = Adapt.adapt_structure(to, sampler.values)
    SamplerWithValues(sampler.total,cumulative_weights,values)
end

function (sampler::SamplerWithValues)(x)
    i = searchsortedfirst(sampler.cumulative_weights, x)
    i = clamp(i, firstindex(sampler.cumulative_weights), lastindex(sampler.cumulative_weights))
    @inbounds sampler.values[i]
end

# TODO: use Alias sampler
struct Sampler{T,A<:AbstractVector{T}}
    total::T
    cumulative_weights::A
end

function Sampler(weights)
    total = sum(weights)
    cumulative_weights = cumsum(weights)
    return Sampler(total,cumulative_weights)
end

function Adapt.adapt_structure(to, sampler::Sampler)
    cumulative_weights = Adapt.adapt_structure(to, sampler.cumulative_weights)
    Sampler(sampler.total,cumulative_weights)
end

function (sampler::Sampler)(x)
    i = searchsortedfirst(sampler.cumulative_weights, x)
    clamp(i, firstindex(sampler.cumulative_weights), lastindex(sampler.cumulative_weights))
end


function sample(rng::AbstractRNG, weights::CuVector,nshots::Integer)
    sampler = Sampler(weights)

    dices = sampler.total .* rand(rng, nshots)

    sampler.(dices)
end

function sample(rng::AbstractRNG, weights::CuVector,nshots::Integer, obs::CuVector)
    sampler = SamplerWithValues(weights,obs)

    dices = sampler.total .* rand(rng, nshots)

    sampler.(dices)
end


N = 20
ntot = 2^(N-4)
space = Subspace(N, sort(randperm(1 << N)[1:ntot] .- 1))
r = SubspaceArrayReg(randn(ComplexF64, ntot), space)
normalize!(r)
dr = cu(r)


weights = abs2.(relaxedvec(dr))

using BenchmarkTools 

@benchmark CUDA.@sync sample(CURAND.default_rng(),weights,10000)

