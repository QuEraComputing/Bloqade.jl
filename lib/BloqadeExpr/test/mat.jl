using Test
using Random
using YaoBlocks
using BloqadeExpr
using Yao
using SparseArrays

@testset "mat" begin
    ss = Subspace(10, sort(randperm(1023)[1:100]))
    U = zeros(ComplexF64, 1024, 1024)
    U[ss, ss] = rand_unitary(length(ss))
    H = U * U'
    for g in [
        put(10, 2 => X),
        put(10, 3 => Rx(0.4)),
        put(10, 2:6 => matblock(rand_unitary(32))) * 2,
        control(10, (3, -5), (2, 7) => matblock(rand_unitary(4))),
        chain(matblock(U), matblock(U)),
        igate(10),
        subroutine(10, kron(X, X), (5, 2)),
        control(10, (3, -5), (2, 7) => matblock(rand_unitary(4)))',
        put(10, 2 => X) + matblock(rand_unitary(1024)),
        time_evolve(matblock(H), 0.3),
    ]
        M = SparseMatrixCSC(mat(g))[ss, ss]
        @test M ≈ mat(g, ss)
    end
end

@testset "XPhase = PdPhase + PuPhase" begin
    @test mat(XPhase(1.0)) ≈ mat(PdPhase(1.0)) + mat(PuPhase(1.0))
end

@testset "Int128" begin
    h = rydberg_h([(rand(2) .* 10...,) for i=1:66], 0.2, 0.3, 0.3, 0.6)
    s = Subspace(127, rand(Int128, 100) .>> 1)
    @test mat(h, s) isa SparseMatrixCSC
end

@testset "3-level operators" begin
    ids_r = [5, 6, 8, 9]
    ids_hf = [1, 2, 4, 5]

    SX_r = mat(SumOfX_1r(2))
    @test isapprox(SX_r[ids_r, ids_r], mat(SumOfX(2)))
    SX_hf = mat(SumOfX_01(2))
    @test SX_hf[ids_hf, ids_hf] ≈ mat(SumOfX(2))
    SZ_r = mat(SumOfZ_1r(2))
    @test SZ_r[ids_r, ids_r] ≈ mat(SumOfZ(2))
    SZ_hf = mat(SumOfZ_01(2))
    @test SZ_hf[ids_hf, ids_hf] ≈ mat(SumOfZ(2))
    SN_r = mat(SumOfN_r(2))
    @test SN_r[ids_r, ids_r] ≈ mat(SumOfN(2))
    SN_hf = mat(SumOfN_1(2))
    @test SN_hf[ids_hf, ids_hf] ≈ mat(SumOfN(2))

    SXPhase_r = mat(SumOfXPhase_1r(2, 1, [0, 0]))
    @test SX_r ≈ SXPhase_r
    
    ϕ = rand()*2π
    @test mat(PuPhase_1r(ϕ)) + mat(PdPhase_1r(ϕ)) ≈ mat(XPhase_1r(ϕ))
    @test mat(PuPhase_01(ϕ)) + mat(PdPhase_01(ϕ)) ≈ mat(XPhase_01(ϕ))
    @test sum(mat.([Pu_1r, Pd_1r, Pu_01, Pd_01])) ≈ [0 1 0; 1 0 1; 0 1 0]
end
