export measure, measure_mis, measure!, measure_mis!

function Yao.measure(
    ::Yao.ComputationalBasis,
    reg::RydbergReg{N, 1},
    ::Yao.AllLocs;
    nshots::Int = 1,
    rng::AbstractRNG = Random.GLOBAL_RNG,
    ) where N
    BitStr64{N}.(sample(rng, vec(reg.subspace), Weights(abs2.(Yao.relaxedvec(reg))), nshots))
end

function Yao.measure!(
    ::Yao.NoPostProcess,
    ::Yao.ComputationalBasis,
    reg::RydbergReg{N, 1},
    ::Yao.AllLocs;
    rng::AbstractRNG = Random.GLOBAL_RNG,
    ) where N
    ind = sample(rng, 1:length(reg.subspace), Weights(abs2.(Yao.relaxedvec(reg))))
    reg.state .= 0
    reg.state[ind] = 1
    return BitStr64{N}(vec(reg.subspace)[ind])
end

function Yao.measure(; nshots=1)
    reg -> Yao.measure(reg; nshots=nshots)
end
