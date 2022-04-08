using Test
using Adapt
using LinearAlgebra
using ComplexArrays
using ComplexArrays: unwarp_complex_array_type

@testset "ComplexArray(...)" begin
    X = ComplexArray(rand(Float64, 10, 2))
    for i in 1:10
        @test X[i] == X.storage[i, 1] + im * X.storage[i, 2]
    end

    X = ComplexArray(rand(Float64, 10, 10, 2))
    for i in 1:10, j in 1:10
        @test X[i, j] == X.storage[i, j, 1] + im * X.storage[i, j, 2]
    end

    X = ComplexArray{Float64, 1}(undef, 10)
    @test size(X) == (10, )
    X = ComplexArray{Float64, 2}(undef, 10, 10)
    @test size(X) == (10, 10)
    X = ComplexArray{Float64, 3}(undef, 10, 10, 10)
    @test size(X) == (10, 10, 10)

    X = ComplexArray{Float64}(undef, 10)
    @test size(X) == (10, )
    X = ComplexArray{Float64}(undef, 10, 10)
    @test size(X) == (10, 10)
    X = ComplexArray{Float64}(undef, 10, 10, 10)
    @test size(X) == (10, 10, 10)
end

@testset "copyto!(::ComplexArray, ::Array)" begin
    X = rand(ComplexF64, 10)
    Y = ComplexVector{Float64}(undef, (10, ))
    copyto!(Y, X)
    @test Y ≈ X

    X = rand(ComplexF64, 10, 10)
    Y = ComplexMatrix{Float64}(undef, (10, 10))
    copyto!(Y, X)
    @test Y ≈ X
end

@testset "similar" begin
    X = ComplexArray{Float64}(undef, 10)
    Y = similar(X, (10, 10))
    @test size(Y) == (10, 10)
    @test eltype(Y) === ComplexF64
    @test Y isa ComplexArray
end

@testset "adapt" begin
    X = rand(ComplexF64, 10, 5)
    Y = adapt(ComplexArray, X)
    @test Y isa ComplexMatrix
    @test X ≈ Y
end

@testset "linalg" begin
    M = rand(10, 10)
    X = adapt(ComplexArray, rand(ComplexF64, 10))
    Y = similar(X)
    mul!(Y, M, X, true, false)
    @test Y ≈ M * X
end

@testset "broadcast" begin
    X = adapt(ComplexArray, rand(ComplexF64, 10))
    @test X .+ X ≈ 2 * X
end

@testset "unwarp_complex_array_type" begin
    X = adapt(ComplexArray, rand(ComplexF64, 10))
    @test unwarp_complex_array_type(X, X, 2.0) isa Tuple{Matrix{Float64}, Matrix{Float64}, Float64}
    @test_throws ErrorException unwarp_complex_array_type(X, rand(2, 2), X)
    @test_throws ErrorException unwarp_complex_array_type(rand(2, 2), X, X)
end
