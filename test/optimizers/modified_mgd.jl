using RydbergEmulator
using RydbergEmulator: modified_mgd, _optimal, multivariate_quadratic
using Test, Random

@testset "test optimization" begin
    Random.seed!(6)
    params = [1.0, 2.0, -2.0, 1.0, 0.0, 1.0]
    _optimal(multivariate_quadratic, params) == [-1.0, 1.0]
end

@testset "test optimization" begin
    Random.seed!(5)
    #rosen(x) = (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2
    #rosen(x) = (1.0 - x[1])^2 + 100.0 * (x[2] - x[1])^2
    rosen(x) = sin(2*x[1]+0.5) + cos(x[2]-0.5) * 0.7 + 0.5*sin(x[1])*cos(x[2])
    x0 = randn(2)
    x = modified_mgd(rosen, x0; δ=0.2, k=50,
                                ϵ=1e-8, n=1000,
                                p0=randn(6), lr=0.5)
    allpass = true
    for i=1:100
        xi = x .+ 0.2 .* randn(2)
        allpass = allpass && rosen(x) < rosen(xi)
    end
    @test allpass
end
