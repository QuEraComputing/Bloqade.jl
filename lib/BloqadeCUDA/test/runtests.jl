using Test

@testset "krylov" begin
    include("krylov.jl")
end

@testset "ode" begin
    include("ode.jl")
end
