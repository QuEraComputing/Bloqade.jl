using RydbergEmulator
using Test
using Random
Random.seed!(42)

include("utils.jl")

@testset "atoms" begin
    include("atoms.jl")
end

@testset "hamiltonian" begin
    include("hamiltonian.jl")
end

@testset "QAOA emulator" begin
    include("emulate.jl")
end

@testset "Yao interfaces" begin
    include("register.jl")
    include("measure.jl")
end

@testset "unit disk" begin
    include("unit_disk.jl")
end
