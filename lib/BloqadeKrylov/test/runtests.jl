using Test
using BloqadeKrylov

if "docstring" in ARGS
    # include("docstrings.jl")
    exit()
end

@testset "expmv" begin
    include("expmv.jl")
end

@testset "emulate" begin
    include("emulate.jl")
end

# @testset "forwarddiff" begin
#     include("forwarddiff.jl")
# end
