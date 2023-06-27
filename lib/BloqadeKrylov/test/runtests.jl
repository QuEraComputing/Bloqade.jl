using Test
using BloqadeKrylov


if "docstring" in ARGS
    # include("docstrings.jl")
    exit()
end

@testset "onenormest" begin
    include("onenormest.jl")

end

#=
@testset "expm_multiply" begin
    include("expm_multiply.jl")
end


@testset "ValHamiltonian" begin
    include("utils.jl")
end
=#

#=
@testset "expmv" begin
    include("expmv.jl")
end

@testset "emulate" begin
    include("emulate.jl")
end


@testset "cfet42_sinX" begin
    include("cfet42_sinX.jl")
end

@testset "cfet65_sinX" begin
    include("cfet65_sinX.jl")
end

@testset "cfet811_sinX" begin
    include("cfet811_sinX.jl")
end


@testset "krylov_sinX" begin
    include("krylov_sinX.jl")
end

@testset "magnus4_sinX" begin
    include("magnus4_sinX.jl")
end


@testset "forwarddiff" begin
    include("forwarddiff.jl")
end
=#