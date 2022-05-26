function YaoAPI.measure(
    ::YaoAPI.ComputationalBasis,
    reg::SubspaceArrayReg,
    ::YaoAPI.AllLocs;
    nshots::Int = 1,
    rng::AbstractRNG = Random.GLOBAL_RNG,
)
    return BitStr64{reg.natoms}.(sample(rng, vec(reg.subspace), Weights(abs2.(relaxedvec(reg))), nshots))
end

function YaoAPI.measure!(
    ::YaoAPI.NoPostProcess,
    ::YaoAPI.ComputationalBasis,
    reg::SubspaceArrayReg,
    ::YaoAPI.AllLocs;
    rng::AbstractRNG = Random.GLOBAL_RNG,
)
    ind = sample(rng, 1:length(reg.subspace), Weights(abs2.(relaxedvec(reg))))
    fill!(reg.state, 0)
    reg.state[ind] = 1
    return BitStr64{reg.natoms}(vec(reg.subspace)[ind])
end

function YaoAPI.measure(; nshots = 1)
    return reg -> YaoAPI.measure(reg; nshots = nshots)
end

# TODO: remove this after https://github.com/QuantumBFS/Yao.jl/issues/338
function YaoAPI.expect(op::AbstractBlock, reg::SubspaceArrayReg)
    return reg' * apply!(copy(reg), op)
end
