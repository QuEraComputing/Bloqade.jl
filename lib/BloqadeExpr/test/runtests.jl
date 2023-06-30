using Test
using BloqadeExpr


if "docstring" in ARGS
    # include("docstrings.jl")
    exit()
end

@testset "linear algebra on linear map: mul!()" begin
    include("linalg_mul.jl")
end

@testset "linear algebra on trace: tr()" begin
    include("linalg_tr.jl")
end

@testset "linear algebra add_I()" begin
    include("linalg_addI.jl")
end

@testset "assertions" begin
    include("assert.jl")
end

@testset "interfaces" begin
    include("interface.jl")
end

@testset "types" begin
    include("types.jl")
end

@testset "lower" begin
    include("lower.jl")
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

@testset "printings" begin
    include("printings.jl")
end

@testset "3-level supports" begin
    include("3-level_supports.jl")
end
