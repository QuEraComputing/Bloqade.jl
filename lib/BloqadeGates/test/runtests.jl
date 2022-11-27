using BloqadeGates
using Test

if "docstring" in ARGS
    # include("docstrings.jl")
    exit()
end

@testset "utils.jl" begin
    include("utils.jl")
end

@testset "types.jl" begin
    include("types.jl")
end

@testset "pulse_generator.jl" begin
    include("pulse_generator.jl")
end

@testset "Predefined Pulses" begin
    include("predefined_pulses.jl")
end

@testset "Pulse Sequence Zoo" begin
    include("predefined_sequences.jl")
end
