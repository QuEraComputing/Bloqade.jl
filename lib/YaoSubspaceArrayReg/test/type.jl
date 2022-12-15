using Test
using Adapt
using Random
using YaoArrayRegister: ArrayReg, probs
using YaoSubspaceArrayReg


@testset "zero_state" begin
    subspace_v = [0, 1, 4, 8]
    raw_st = zeros(ComplexF64, length(subspace_v))
    raw_st[1] = 1
    r = zero_state(Subspace(10, subspace_v))
    @test state(r) ≈ raw_st
    @test YaoSubspaceArrayReg.basis(r) == YaoSubspaceArrayReg.BitStr{10}.(subspace_v)

    raw_st = rand(ComplexF64, length(subspace_v))
    r = SubspaceArrayReg(raw_st, Subspace(10, subspace_v))

    @test adapt(Array{ComplexF32}, r).state isa Array{ComplexF32}

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
    space = Subspace(5, sort(randperm(1 << 5)[1:6] .- 1))
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
    a = SubspaceArrayReg(ComplexF64[0.0, 0.8, 0.6, 0.0], s)
    b = rand_state(s)
    @test a == a
    @test most_probable(a, 2) == YaoSubspaceArrayReg.BitStr{5}.([1, 4])
    @test YaoSubspaceArrayReg.regadd!(copy(a), b).state ≈ (a.state + b.state)
    @test YaoSubspaceArrayReg.regsub!(copy(a), b).state ≈ (a.state - b.state)
    @test YaoSubspaceArrayReg.regscale!(copy(a), 0.5).state ≈ 0.5 * a.state
    @test (a * 0.5).state ≈ 0.5 * a.state
    @test (-a).state ≈ -a.state
    @test (a / 0.5).state ≈ a.state / 0.5
    @test (a + b).state ≈ a.state + b.state
    @test (a - b).state ≈ a.state - b.state
end

@testset "ArrayReg(::SubspaceArrayReg)" begin
    subspace_v = [0, 1, 4, 8]
    space = Subspace(10, subspace_v)
    r = rand_state(space)
    @test ArrayReg(r).state[subspace_v.+1] ≈ r.state
end

@testset "probs(::SubspaceArrayReg)" begin
    subspace_v = [0, 1, 4, 8]
    space = Subspace(10, subspace_v)
    @test sum(probs(rand_state(space))) ≈ 1.0
end

@testset "getindex(::SubspaceArrayReg, ::DisStr)" begin
    subspace = [0, 1, 4, 8]
    s = Subspace(5, subspace)
    reg = zero_state(s)
    @test reg[bit"00000"] == 1
    @test reg[bit"00001"] == 0
    @test reg[bit"00011"] == 0
    @test reg[bit"00100"] == 0
    @test reg[bit"01000"] == 0
    @test reg[bit"01001"] == 0
    @test reg[bit"01011"] == 0
    @test reg[bit"01111"] == 0
    @test reg[bit"11111"] == 0
end
