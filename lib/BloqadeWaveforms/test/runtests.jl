using Test

@testset "waveform" begin
    include("waveform.jl")
end

@testset "smooth" begin
    include("smooth.jl")
end

@testset "descretize" begin
    include("discretize.jl")
end