using Test
using CUDA
using Adapt
using EaRydCUDA
using ComplexArrays
using LinearAlgebra

@testset "rmul(::ComplexArray{.., <:CuArray}, ::Complex)" begin
    X = rand(ComplexF64, 10)
    Y = adapt(ComplexArray, X)
    dY = cu(Y)

    alpha = 1.3+2.1im
    rmul!(dY, alpha)
    @test adapt(Array, dY) ≈ X * alpha

    # test on multiple blocks
    X = rand(ComplexF64, 10000)
    Y = adapt(ComplexArray, X)
    dY = cu(Y)
    alpha = 1.3+2.1im
    rmul!(dY, alpha)

    @test adapt(Array, dY) ≈ X * alpha
end

alpha = 0.09359693561994409 - 0.0im
X = rand(ComplexF64, 1000)
Y = adapt(ComplexArray, X)
dY = adapt(CuArray, Y)
rmul!(dY, alpha)
adapt(Array, dY).storage ≈ (Y * alpha).storage
Y * alpha
