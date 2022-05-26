using Test
using InteractiveUtils
using BloqadeExpr: is_time_function, assert_nsites

@testset "is_time_function" begin
    @test is_time_function(sin) == true
    @test is_time_function(push!) == false
    tf_1(::Float64) = 0.0
    @test is_time_function(tf_1) == true
    tf_2(::Complex) = 0.0
    @test is_time_function(tf_2) == false
end

@testset "assert_nsites" begin
    @test assert_nsites(5, [1, 2, 3, 4, 5], :test) === nothing
    @test assert_nsites(5, 2.0, :test) === nothing
    @test_throws ArgumentError assert_nsites(4, [1, 2, 3, 4, 5], :test)
end
