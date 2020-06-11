export measure, measure_mis, measure!, measure_mis!

function Yao.measure(
    ::Yao.ComputationalBasis,
    reg::RydbergReg{N, 1},
    ::Yao.AllLocs;
    nshots::Int = 1,
    rng::AbstractRNG = Random.GLOBAL_RNG,
    ) where N
    sample(rng, reg.subspace, Weights(abs2.(relaxedvec(reg))), nshots)
end

function Yao.measure!(
    ::Yao.NoPostProcess,
    ::Yao.ComputationalBasis,
    reg::RydbergReg{N, 1},
    ::Yao.AllLocs;
    rng::AbstractRNG = Random.GLOBAL_RNG,
    ) where N
    ind = sample(rng, 1:length(reg.subspace), Weights(abs2.(relaxedvec(reg))))
    reg.state .= 0
    reg.state[ind] = 1
    return reg.subspace[ind]
end

function Yao.measure(; nshots=1)
    reg -> measure(reg; nshots=nshots)
end
