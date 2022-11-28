using BloqadeGates
using BloqadeExpr
using BloqadeODE, BloqadeKrylov
using Yao
using Test

using BloqadeGates: two_level_indices

@testset "Local pulses" begin
    @testset "C_{$n}Z" for n = 1:5
        atoms = [(cos(2π*i/n), sin(2π*i/n)) for i = 1:n]
        if n == 1
            @test_throws AssertionError local_CkZ(atoms, [i for i = 2:n], [1])
        else
            seq = local_CkZ(atoms, [i for i = 2:n], [1])
            m = mat(seq)
            ids = two_level_indices(n)
            m = m[ids, ids]
            @test operator_fidelity(matblock(m), control(n, collect(2:n), 1=>Z)) ≈ 1
        end
    end
    
    @testset "C_{$n}NOT" for n = 1:5
        atoms = [(cos(2π*i/n), sin(2π*i/n)) for i = 1:n]
        if n == 1
            @test_throws AssertionError local_CkZ(atoms, [i for i = 2:n], [1])
        else
            seq = local_CkNOT(atoms, [i for i = 2:n], [1])
            m = mat(seq)
            ids = two_level_indices(n)
            m = m[ids, ids]
            @test operator_fidelity(matblock(m), control(n, collect(2:n), 1=>X)) ≈ 1
        end
    end
end

@testset "The Levine-Pichler gate" begin
    atoms = [(0.0, 0.0), (0.0, 4.0)]
    seq = global_levine_pichler(atoms)
    push!(seq, global_single_qubit_gate(atoms, Rz(2π - 2.38076)))
    m = mat(seq)
    ids = two_level_indices(2)
    m = m[ids, ids]
    @test 1-1e-6 < operator_fidelity(matblock(m), control(2, 1, 2=>Z)) < 1+1e-6
end
