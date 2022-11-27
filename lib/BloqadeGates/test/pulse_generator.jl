using BloqadeGates
using Test
using BloqadeGates: apply_mask

@testset "apply_mask" begin
    N = 3
    mask = rand(N)
    fs = apply_mask(sin, mask)
    for i = 1:N
        @test (fs[i])(1) ≈ mask[i]*sin(1)
    end
    @test_throws ErrorException apply_mask([1, 2, 3], mask)
end