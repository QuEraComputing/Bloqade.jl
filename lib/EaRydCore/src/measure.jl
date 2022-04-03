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

# TODO: remove this after https://github.com/QuantumBFS/Yao.jl/issues/338
function Yao.expect(op::AbstractBlock, reg::RydbergReg)
    return reg' * apply!(copy(reg), op)
end

function YaoBlocks.regadd!(reg::RydbergReg, reg2::RydbergReg)
    @assert reg.natoms == reg2.natoms
    @assert reg.layout == reg2.layout
    @assert length(reg.subspace) == length(reg2.subspace)
    reg.state .+= reg2.state
    return reg
end
