using RydbergEmulator
using Test
using Random
Random.seed!(42)

@testset "Yao interfaces" begin
    include("register.jl")
    include("measure.jl")
end

@testset "QAOA emulator" begin
    include("qaoa.jl")
end

@testset "qaoa_mis" begin
    include("qaoa_mis.jl")
end

@testset "RydbergEmulator.jl" begin
    include("unit_disk.jl")
    include("hamiltonian.jl")
end
