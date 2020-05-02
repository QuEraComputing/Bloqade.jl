using LinearAlgebra
using Test
using RydbergEmulator

@testset "CMA-ES RosenBrock" begin
	# Testing: CMA-ES
	N = 2
    rosenbrock(x) = (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2
	# termination condition
	result = cmaes(rosenbrock, randn(2); num_offsprings=12, num_parents=3, maxiter=10000, tol=1e-10)
    @test length(result) == N
	@test ≈(result, ones(N), atol=1e-5)
    @test ≈(rosenbrock(result), 0.0, atol=1e-5)
end
