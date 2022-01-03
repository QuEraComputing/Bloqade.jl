using Test
using EaRydCore
using Yao

@testset "zero_state" begin
    subspace_v = [0, 1, 4, 8]
    raw_st = zeros(ComplexF64, length(subspace_v))
    raw_st[1] = 1
    @test state(zero_state(Subspace(10, subspace_v))) ≈ raw_st

    raw_st = rand(ComplexF64, length(subspace_v))
    r = RydbergReg(raw_st, Subspace(10, subspace_v))
    set_zero_state!(r)
    @test r.state[1] == 1
    @test relaxedvec(r) isa Vector

    raw_st = rand(ComplexF64, length(subspace_v))
    r = RydbergReg(raw_st, Subspace(10, subspace_v))
    set_zero_state!(r)
    @test r.state[1] == 1

    @test nactive(r) == 10
    @test nqubits(r) == 10
    @test relaxedvec(r) isa Vector
    @test state(r) isa Vector
    @test statevec(r) isa Vector
    @test isnormalized(r)

    @test set_zero_state!(Yao.rand_state(5)) ≈ Yao.zero_state(5)
    @test_throws DimensionMismatch RydbergReg(rand(5), Subspace(5, rand(Int, 6)))

    # RealLayout
    r = zero_state(Subspace(5, subspace_v), RealLayout())
    @test r.state[1, 1] == 1
    @test r.state[1, 2] == 0
    @test all(r.state[2:end, :] .== 0)
end

@testset "rand_state" begin
    subspace = [0, 1, 4, 8]
    s = Subspace(5, subspace)
    @test isnormalized(rand_state(s))

    @test isnormalized(rand_state(s, RealLayout()))
end

@testset "product_state" begin
    subspace = [0, 1, 4, 8]
    s = Subspace(5, subspace)
    @test_throws ErrorException product_state(bit"010", s)
    r = product_state(bit"1000", s)
    @test isone(r.state[end])
end

@testset "constructors" begin
    r = rand_state(Subspace(10, [0, 1, 4, 8]), ComplexLayout())
    @test RydbergReg{RealLayout}(r) isa RydbergReg{RealLayout}
    @test RydbergReg{ComplexLayout}(r) isa RydbergReg{ComplexLayout}

    r = rand_state(Subspace(10, [0, 1, 4, 8]), RealLayout())
    @test RydbergReg{RealLayout}(r) isa RydbergReg{RealLayout}
    @test RydbergReg{ComplexLayout}(r) isa RydbergReg{ComplexLayout}
end
