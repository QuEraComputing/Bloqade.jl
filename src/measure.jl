export measure, measure_mis, measure!, measure_mis!, expect_mis

function YaoBase.measure(
    ::ComputationalBasis,
    reg::RydbergReg{1},
    ::AllLocs;
    nshots::Int = 1,
    rng::AbstractRNG = Random.GLOBAL_RNG,
    )
    sample(rng, reg.subspace, Weights(abs2.(relaxedvec(reg))), nshots)
end

function YaoBase.measure!(
    ::YaoBase.NoPostProcess,
    ::ComputationalBasis,
    reg::RydbergReg{1},
    ::AllLocs;
    rng::AbstractRNG = Random.GLOBAL_RNG,
    )
    ind = sample(rng, 1:length(reg.subspace), Weights(relaxedvec(reg)))
    reg.state .= 0
    reg.state[ind] = 1
    return reg.subspace[ind]
end

"""
    measure_mis!(reg::RydbergReg)

Measure the independent set size on `reg`, and collapse the state.
"""
function measure_mis!(reg::RydbergReg)
    count_ones(measure!(reg))
end

"""
    measure_mis(reg::RydbergReg; nshots=1)

Measure the independent set size on `reg` for `nshots` times.
"""
function measure_mis(reg::RydbergReg; nshots=1)
    count_ones.(measure(reg; nshots=nshots))
end

"""
    expect_mis(reg::RydbergReg)

Get the expected size of independent set on `reg`.
"""
function expect_mis(reg::RydbergReg)
    sum(t -> abs2(t[2]) * count_ones(t[1]), zip(reg.subspace, relaxedvec(reg)))
end
