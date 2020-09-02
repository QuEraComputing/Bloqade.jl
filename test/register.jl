using Test
using RydbergEmulator
using Yao

@testset "rydberg register" begin
    subspace = [0, 1, 4, 8]
    raw_st = zeros(ComplexF64, length(subspace), 5)
    raw_st[1, :] .= 1
    @test state(RydbergEmulator.zero_state(10, Subspace(subspace); nbatch=5)) ≈ raw_st

    raw_st = rand(ComplexF64, length(subspace), 5)
    r = RydbergReg{10}(raw_st, Subspace(subspace))
    set_zero_state!(r)
    @test all(r.state[1, :] .== 1)
    @test relaxedvec(r) isa Matrix

    raw_st = rand(ComplexF64, length(subspace))
    r = RydbergReg{10}(raw_st, Subspace(subspace))
    set_zero_state!(r)
    @test all(r.state[1, :] .== 1)

    @test nactive(r) == 10
    @test nqubits(r) == 10
    @test relaxedvec(r) isa Vector
    @test state(r) isa Matrix
    @test statevec(r) isa Vector
    @test isnormalized(r)

    @test set_zero_state!(Yao.rand_state(5)) ≈ Yao.zero_state(5)

    @test_throws DimensionMismatch RydbergReg{5}(rand(5), Subspace(rand(Int, 6)))
end

@testset "rand_state" begin
    subspace = [0, 1, 4, 8]
    s = Subspace(subspace)
    @test isnormalized(rand_state(5, s)) 
end

@testset "product_state" begin
    subspace = [0, 1, 4, 8]
    s = Subspace(subspace)
    @test_throws ErrorException product_state(5, bit"010", s)
    r = product_state(5, bit"1000", s)
    @test isone(r.state[end])    
end
