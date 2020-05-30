export measure, measure_mis, measure!, measure_mis!

function YaoBase.measure(
    ::ComputationalBasis,
    reg::RydbergReg{N, 1},
    ::AllLocs;
    nshots::Int = 1,
    rng::AbstractRNG = Random.GLOBAL_RNG,
    ) where N
    sample(rng, reg.subspace, Weights(abs2.(relaxedvec(reg))), nshots)
end

function YaoBase.measure!(
    ::YaoBase.NoPostProcess,
    ::ComputationalBasis,
    reg::RydbergReg{N, 1},
    ::AllLocs;
    rng::AbstractRNG = Random.GLOBAL_RNG,
    ) where N
    ind = sample(rng, 1:length(reg.subspace), Weights(abs2.(relaxedvec(reg))))
    reg.state .= 0
    reg.state[ind] = 1
    return reg.subspace[ind]
end

function YaoBase.measure(; nshots=1)
    reg -> measure(reg; nshots=nshots)
end
