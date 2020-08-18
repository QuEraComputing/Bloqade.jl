using Test
using RydbergEmulator
using Yao

@testset "rydberg register" begin
    subspace = [0, 1, 4, 8]
    raw_st = zeros(ComplexF64, length(subspace), 5)
    raw_st[1, :] .= 1
    @test state(RydbergEmulator.zero_state(10, Subspace(subspace); nbatch=5)) ≈ raw_st

    raw_st = rand(ComplexF64, length(subspace), 5)
    r = RydbergEmulator.RydbergReg{10}(raw_st, Subspace(subspace))
    set_zero_state!(r)
    @test all(r.state[1, :] .== 1)
end
