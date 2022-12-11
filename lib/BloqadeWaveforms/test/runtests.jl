using Test

if "docstring" in ARGS
    include("docstrings.jl")
    exit()
end

@testset "waveform" begin
    include("waveform.jl")
end

@testset "smooth" begin
    include("smooth.jl")
end

@testset "interpolate" begin
    include("interpolate.jl")
end