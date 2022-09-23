using Test
using BloqadeExpr, YaoBlocks, YaoAPI, BitBasis
using SparseArrays
using LinearAlgebra

@testset "getindex" begin
    pb = rydberg_h([(5*randn(), 5*randn()) for i=1:5]; Ω=0.3, Δ=0.5)
    mpb = mat(pb)
    allpass = true
    for i=basis(pb), j=basis(pb)
        allpass &= pb[i, j] ≈ mpb[Int(i)+1, Int(j)+1]
    end
    @test allpass

    allpass = true
    for j=basis(pb)
        allpass &= vec(pb[:, j]) ≈ mpb[:, Int(j)+1]
        allpass &= vec(pb[j,:]) ≈ mpb[Int(j)+1,:]
        allpass &= vec(pb[:, EntryTable([j], [1.0+0im])]) ≈ mpb[:, Int(j)+1]
        allpass &= vec(pb[EntryTable([j], [1.0+0im]),:]) ≈ mpb[Int(j)+1,:]
    end
    @test allpass
end

@testset "isreal" begin

    @test isreal(SumOfN(nsites=2))

    @test isreal(SumOfZ(nsites=2))

    scale = 2 * SumOfX(nsites=3)
    @test isreal(scale)
    scale = im * SumOfX(nsites=3)
    @test !isreal(scale)
end

@testset "Hamiltonian Term Default Values" begin
    pauli_zs_term = SumOfZ(2)
    @test pauli_zs_term.Δ == 1
    
    number_ops_term = SumOfN(2)
    @test number_ops_term.Δ == 1

    pauli_xs_term = SumOfX(2)
    @test pauli_xs_term.Ω == 1
end

@testset "Qudit Numbers" begin
    @test YaoAPI.nqudits(PdPhase(1.0)) == 1
    @test YaoAPI.nqudits(PuPhase(1.0)) == 1
    @test YaoAPI.nqudits(SumOfZ(10)) == 10
end

@testset "SumOfZ Equality Operator" begin
    @test SumOfZ(10, 1.0) == SumOfZ(10, 1.0)
    @test !(SumOfZ(1, 1.0) == SumOfZ(2, 1.0))
    @test !(SumOfZ(5, 0.5) == SumOfZ(5, 0.3))
end

@testset "Hamiltonian Size" begin
    hamiltonian = BloqadeExpr.Hamiltonian(Float64, SumOfX(6, sin) + SumOfZ(6, cos))
    @test size(hamiltonian) == (64, 64)
    @test size(hamiltonian, 2) == 64
    
    step_hamiltonian = hamiltonian(0.1)
    @test size(step_hamiltonian) == (64, 64)
    @test size(step_hamiltonian, 2) == 64
end

@testset "Step Hamiltonian Norm" begin 
    hamiltonian = BloqadeExpr.Hamiltonian(Float64, SumOfZ(1, sin))
    step_hamiltonian = hamiltonian(π/2)
    @test LinearAlgebra.opnorm(step_hamiltonian) == 1.0
end
