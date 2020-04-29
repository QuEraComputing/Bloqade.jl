using RydbergEmulator
using Test
using Random
Random.seed!(42)

@testset "measure" begin
    include("measure.jl")
end

@testset "optimizers" begin
    include("optimizers/optimizers.jl")
end

@testset "qaoaopt" begin
    include("qaoaopt.jl")
end

@testset "RydbergEmulator.jl" begin
    include("unit_disk.jl")
    include("hamiltonian.jl")
end
