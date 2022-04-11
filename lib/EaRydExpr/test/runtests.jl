using Test
using EaRydExpr

@testset "assertions" begin
    include("assert.jl") 
end

@testset "interfaces" begin
    include("interface.jl")
end

@testset "types" begin
    include("types.jl")
end

@testset "linear algebra on linear map" begin
    include("linalg.jl")    
end

@testset "matrix construction" begin
    include("mat.jl")    
end

@testset "space" begin
    include("space.jl")    
end

@testset "units" begin
    include("units.jl")
end

@testset "atoms" begin
    include("atoms.jl")
end
