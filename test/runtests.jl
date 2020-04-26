using RydbergEmulator
using Test
using Random
Random.seed!(42)

@testset "optimizers" begin
    include("optimizers/optimizers.jl")
end

@testset "RydbergEmulator.jl" begin
    include("unit_disk.jl")
    include("hamiltonian.jl")
end
