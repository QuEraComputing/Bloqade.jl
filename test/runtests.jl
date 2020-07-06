using RydbergEmulator
using Test
using Random
Random.seed!(42)

include("utils.jl")

@testset "hamiltonian" begin
    include("hamiltonian.jl")
end

@testset "QAOA emulator" begin
    include("qaoa.jl")
end

@testset "Yao interfaces" begin
    include("register.jl")
    include("measure.jl")
end

# TODO: move this to StochasitcOptimizer
# @testset "qaoa_mis" begin
#     include("qaoa_mis.jl")
# end

@testset "unit disk" begin
    include("unit_disk.jl")
end
