export measure, measure_mis, measure!, measure_mis!

function Yao.measure(
    ::Yao.ComputationalBasis,
    reg::RydbergReg{ComplexLayout},
    ::Yao.AllLocs;
    nshots::Int = 1,
    rng::AbstractRNG = Random.GLOBAL_RNG,
    )
    BitStr64{reg.natoms}.(
        sample(
            rng,
            vec(reg.subspace),
            Weights(abs2.(Yao.relaxedvec(reg))),
            nshots
        )
    )
end

function Yao.measure!(
    ::Yao.NoPostProcess,
    ::Yao.ComputationalBasis,
    reg::RydbergReg,
    ::Yao.AllLocs;
    rng::AbstractRNG = Random.GLOBAL_RNG,
    )
    ind = sample(rng, 1:length(reg.subspace), Weights(abs2.(Yao.relaxedvec(reg))))
    fill!(reg.state, 0)
    reg.state[ind] = 1
    return BitStr64{reg.natoms}(vec(reg.subspace)[ind])
end

function Yao.measure(; nshots=1)
    reg -> Yao.measure(reg; nshots=nshots)
end
