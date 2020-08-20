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
end
