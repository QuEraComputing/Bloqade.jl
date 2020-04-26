using LsqFit, RydbergEmulator
using StaticArrays, NiLang
using Test, Random

@testset "lsq fit" begin
    @. model(x, p) = p[1] + p[2]*x + p[3]*x^2
    xdata = randn(15) * 10
    ydata = 0.5 .* xdata .^ 2 .+ xdata .* 2.1 .+ 0.3
    p0 = zeros(3)
    fit = curve_fit(model, xdata, ydata, p0)
    @test fit.param ≈ [0.3, 2.1, 0.5]

    @. model(x, p) = p[1] + p[2]*x[1,:] + p[3]*x[2,:] + p[4]*x[1,:]^2 + p[5]*x[2,:]^2 + p[6]*x[1,:]*x[2,:]
    xdata = randn(2, 50) * 10
    q = [0.0, 2.1, 0.3, 0.5, 1.0, 4.0]
    ydata = @. q[1] + q[2]*xdata[1,:] + q[3]*xdata[2,:] + q[4]*xdata[1,:]^2 + q[5]*xdata[2,:]^2 + q[6]*xdata[1,:]*xdata[2,:]
    p0 = zeros(6)
    fit = curve_fit(model, xdata, ydata, p0)
    @test fit.param ≈ q

    @. model(x, p) = p[1] + p[2]*getindex(x, 1) + p[3]*getindex(x,2) + p[4]*getindex(x,1)^2 + p[5]*getindex(x,2)^2 + p[6]*getindex(x,1)*getindex(x,2)
    xdata = randn(MVector{2,Float64}, 50) * 10
    q = [0.0, 2.1, 0.3, 0.5, 1.0, 4.0]
    ydata = @. q[1] + q[2]*getindex(xdata, 1) + q[3]*getindex(xdata,2) + q[4]*getindex(xdata,1)^2 + q[5]*getindex(xdata,2)^2 + q[6]*getindex(xdata,1)*getindex(xdata,2)
    p0 = zeros(6)
    fit = curve_fit(model, xdata, ydata, p0)
    @test fit.param ≈ q
end

@testset "multi_variate_quadratic" begin
    x = randn(2)
    p = randn(6)
    res = RydbergEmulator.multivariate_quadratic(0.0, x, p)[1]
    @test res ≈ p[1] + p[2]*x[1] + p[3]*x[2] + p[4]*x[1]^2 + p[5]*x[1]*x[2] + p[6]*x[2]^2
    _, _, gx, gp = RydbergEmulator.multivariate_quadratic'(Val(1), 0.0, x, p)
    @test NiLang.AD.grad(gx) ≈ [p[2]+p[4]*2*x[1]+p[5]*x[2], p[3]+p[6]*2*x[2]+p[5]*x[1]]
    @test NiLang.AD.grad(gp) ≈ [1, x[1], x[2], x[1]^2, x[1]*x[2], x[2]^2]
end

@testset "test optimization" begin
    Random.seed!(6)
    #rosen(x) = (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2
    rosen(x) = (1.0 - x[1])^2 + 10.0 * (x[2] - x[1])^2
    x0 = [1.0, 1.0] .+ randn(2) .* 0.2
    @show rosen(x0)
    x = mgd(rosen, x0; γ=0.5, δ=0.6, k=10, α=1.0, A=2.0, ξ=1.5,
                                ϵ=1e-8, n=1000,
                                p0=randn(6))
    @test isapprox(rosen(x), 0; atol=1e-3)
end
