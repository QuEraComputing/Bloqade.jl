using Test, BloqadeQMC
using Statistics

if "docstring" in ARGS
    # include("docstrings.jl")
    exit()
end

@testset "Probability Vectors" begin
    @testset "Small Vector of size $n" for n in 1:10
        p = float.(rand(1:10, n))
        mean_ = p' * collect(1:size(p,1)) / sum(p)
        var_ = (p' * (collect(1:size(p,1)) .^ 2) / sum(p)) - mean_^2

        @testset "ProbabilityAlias" begin
            pvec = ProbabilityAlias(p)
            X = [rand(pvec) for _ in 1:1_000_000]

            @test isapprox(mean_, mean(X), atol=0.05)
            @test isapprox(var_, var(X), atol=0.05)
        end
    end
end

H = @inferred LTFIM((10,), 1.0, 1.0, 1.0)
@inferred BinaryGroundState BinaryGroundState(H, 1000)
@inferred BinaryThermalState BinaryThermalState(H, 1000)


THRESHOLD = 2.576  # 99% Two-sided CI of the t-distribution with infinite dofs

include("rydberg.jl")
include("interface.jl")

# include("ltfim.jl")
# include("tfim.jl")
