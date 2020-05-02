using Test, Random
using RydbergEmulator

@testset "SPSA" begin
    rosenbrock(x) = (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2
	Random.seed!(2)
	x0 = randn(2)
	@show rosenbrock(x0)
    result = spsa(rosenbrock, x0, bounds=(-1, 2),
        ac=(0.2, 0.1), maxiter=20000)
    @test length(result) == 2
	@test ≈(result, ones(2), atol=1e-2)
    @test ≈(rosenbrock(result), 0.0, atol=1e-3)
end

@testset "SPSA 2nd order" begin
    rosenbrock(x) = (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2
	Random.seed!(3)
	x0 = randn(2)
	@show rosenbrock(x0)
	# second order function
    result = spsa(rosenbrock, randn(2), bounds=(-1, 2),
        ac=(10.0, 0.03), maxiter=2000, secondorder=true)
    @test length(result) == 2
	@test ≈(result, ones(2), atol=1e-2)
    @test ≈(rosenbrock(result), 0.0, atol=1e-3)
end
