using BloqadeGates
using Yao
using Test

using BloqadeGates: two_level_indices
using Yao.EasyBuild: SqrtX, SqrtY, SqrtW
using Yao.YaoBlocks.ConstGate: I2, X, Y, Z, H, S, T

const single_qubit_gates = [
    I2, phase(rand()*2π - π),
    X, Y, Z, H, S, T, S', T',
    shift(rand()*2π - π), phase(rand()*2π - π), Rx(rand()*2π - π), Ry(rand()*2π - π), Rz(rand()*2π - π),
    SqrtX, SqrtY, SqrtW, 
    matblock(rand_unitary(2))
]

@testset "superposition states 1/(√2^$n) ∑|i⟩" for n = 1:5
    atoms = [(cos(2π*i/n), sin(2π*i/n)) for i = 1:n]
    hs = global_single_qubit_gate(atoms, H)
    reg = zero_state(n; nlevel = 3)
    reg |> hs
    ids = two_level_indices(n)
    @test all(isapprox.(probs(reg)[ids], 1/2^n; atol = 1e-6))
end

@testset "Local single qubit gate" begin
    atoms = [(0.0, 0.0), (4.0, 0.0)]
    n = length(atoms)
    ids = two_level_indices(n)
    for g in single_qubit_gates, loc in 1:n
        p = local_single_qubit_gate(atoms, [loc], g)
        m_p = mat(p)[ids, ids]
        fid = operator_fidelity(matblock(m_p), put(n, loc=>g))
        @test isapprox(fid, 1; atol = 1e-9)
    end
end

@testset "Global single qubit gate" begin
    atoms = [(0.0, 0.0), (4.0, 0.0)]
    n = length(atoms)
    ids = two_level_indices(n)
    for g in single_qubit_gates
        p  = global_single_qubit_gate(atoms, g)
        m_p = mat(p)[ids, ids]
        fid = operator_fidelity(matblock(m_p), kron(g for loc in 1:n))
        @test isapprox(fid, 1; atol = 1e-9)
    end
end