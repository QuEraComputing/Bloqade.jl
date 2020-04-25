using RydbergEmulator
using Test
using Random
Random.seed!(42)

@testset "RydbergEmulator.jl" begin
    include("unit_disk.jl")
    include("hamiltonian.jl")
end
