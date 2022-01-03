using EaRydKrylovEvolution
using Test
using Random
Random.seed!(42)

include("utils.jl")

@testset "expmv" begin
    include("expmv.jl")
end

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
    include("instructs.jl")
    include("mat.jl")
end

@testset "unit disk" begin
    include("unit_disk.jl")
end

@testset "utils" begin
    @test_throws ErrorException EaRydKrylovEvolution.unsafe_log2i(2.2)
end

@testset "serialize" begin
    include("serialize.jl")
end

@testset "mis" begin
    include("mis.jl")
end
