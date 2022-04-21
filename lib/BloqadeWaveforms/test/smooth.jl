using Test
using BloqadeWaveforms
using BloqadeWaveforms.Kernels
using BloqadeWaveforms: edge_pad

@testset "kernels" begin
    @test Kernels._in_radius(0.1, true)
    @test Kernels._in_radius(1.1, true) ≈ 0.0

    @test biweight(0.1) ≈ 0.91884375
    @test biweight(1.1) ≈ 0.0
    @test cosine(0.1) ≈ π/4 * cos(π/2 * 0.1)
    @test gaussian(0.1) ≈ 0.3969525474770118
    @test logistic(0.1) ≈ 0.24937604019289197
    @test parabolic(0.1) ≈ 0.7424999999999999
    @test sigmoid(0.1) ≈ 0.3167249413498117
    @test triangle(0.1) ≈ 0.9
    @test tricube(0.1) ≈ 0.8616075299999999
    @test uniform(0.1) ≈ 1.0
end

@testset "kernel smoother" begin
    wf = piecewise_linear(clocks=[0.0, 2.0, 3.0, 4.0], values=[0.0, 3.0, 1.1, 2.2])
    swf = smooth(wf; kernel_radius=0.5)
    @test abs(swf(0.1) - swf(0.1+1e-3)) ≤ 1e-5        
end
