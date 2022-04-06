using Test
using EaRydKrylov

@testset "expmv" begin
    include("expmv.jl")    
end

@testset "emulate" begin
    include("emulate.jl")    
end

@testset "forwarddiff" begin
    include("forwarddiff.jl")
end
