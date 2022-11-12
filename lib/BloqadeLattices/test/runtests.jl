using Test
using BloqadeLattices

if "docstring" in ARGS
    include("docstrings.jl")
    exit()
end

@testset "lattice" begin
    include("lattice.jl")
end

@testset "region" begin
    include("region.jl")
end

@testset "bounded_lattice" begin
    include("bounded_lattice.jl")
end

@testset "interact" begin
    include("interact.jl")
end

@testset "neighbors" begin
    include("neighbors.jl")
end

@testset "visualize" begin
    include("visualize.jl")
end