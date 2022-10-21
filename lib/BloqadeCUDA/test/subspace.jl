using Test
using Random
using CUDA
using YaoSubspaceArrayReg
using CuYao

function sample(rng::AbstractRNG, wv::AbstractWeights)
    1 == firstindex(wv) ||
        throw(ArgumentError("non 1-based arrays are not supported"))
    t = rand(rng) * sum(wv)
    n = length(wv)
    i = 1
    cw = wv[1]
    while cw < t && i < n
        i += 1
        @inbounds cw += wv[i]
    end
    return i
end

function sample(rng::Union{RNG,GPUArrays.RNG}, weights::CuVector)
    dices = rand(rng, length(weights))
    slots = cumsum(weights)
    return map(dices) do d
        # TODO: use binary search
        findfirst(d, slots)
        searchsortedfirst
    end
end



space = Subspace(10, sort(randperm(1 << 10)[1:76] .- 1))
r = SubspaceArrayReg(randn(ComplexF64, 76), space)

dr = cu(r)
@which measure(r)

@which measure(ComputationalBasis(), r, AllLocs(); nshots=10)

measure(dr)
