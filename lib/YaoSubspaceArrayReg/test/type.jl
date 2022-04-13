using Test
using Random
using YaoSubspaceArrayReg

@testset "zero_state" begin
    subspace_v = [0, 1, 4, 8]
    raw_st = zeros(ComplexF64, length(subspace_v))
    raw_st[1] = 1
    @test state(zero_state(Subspace(10, subspace_v))) ≈ raw_st

    raw_st = rand(ComplexF64, length(subspace_v))
    r = SubspaceArrayReg(raw_st, Subspace(10, subspace_v))
    set_zero_state!(r)
    @test r.state[1] == 1
    @test relaxedvec(r) isa Vector

    raw_st = rand(ComplexF64, length(subspace_v))
    r = SubspaceArrayReg(raw_st, Subspace(10, subspace_v))
    set_zero_state!(r)
    @test r.state[1] == 1

    @test nactive(r) == 10
    @test nqubits(r) == 10
    @test relaxedvec(r) isa Vector
    @test state(r) isa Vector
    @test statevec(r) isa Vector
    @test isnormalized(r)

    @test set_zero_state!(rand_state(5)) ≈ zero_state(5)
    space = Subspace(5, sort(randperm(1<<5)[1:6]))
    @test_throws DimensionMismatch SubspaceArrayReg(rand(5), space)
end

@testset "space" begin
    @test space(rand_state(5)) === fullspace
    @test space(rand_state(Subspace(5, [0, 1, 4, 8]))) isa Subspace
end

@testset "rand_state" begin
    subspace = [0, 1, 4, 8]
    s = Subspace(5, subspace)
    @test isnormalized(rand_state(s))
end

@testset "product_state" begin
    subspace = [0, 1, 4, 8]
    s = Subspace(5, subspace)
    @test_throws ErrorException product_state(bit"010", s)
    r = product_state(bit"1000", s)
    @test isone(r.state[end])

    r = product_state(bit"1000", s)
    @test isone(r' * r)
end

@testset "reg operations" begin
    subspace = [0, 1, 4, 8]
    s = Subspace(5, subspace)
    a = rand_state(s)
    b = rand_state(s)
    @test a == a
    @test YaoSubspaceArrayReg.regadd!(copy(a), b).state ≈ (a.state + b.state)
    @test YaoSubspaceArrayReg.regsub!(copy(a), b).state ≈ (a.state - b.state)
    @test YaoSubspaceArrayReg.regscale!(copy(a), 0.5).state ≈ 0.5 * a.state
    @test (a * 0.5).state ≈ 0.5 * a.state
    @test (-a).state ≈ -a.state
    @test (a / 0.5).state ≈ a.state / 0.5
    @test (a + b).state ≈ a.state + b.state
    @test (a - b).state ≈ a.state - b.state
end