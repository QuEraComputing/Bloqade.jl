using Test, BloqadeQMC
using Statistics
using Random
using RandomNumbers
using Measurements
using BloqadeQMC: bootstrap
using Bootstrap: BasicSampling

rng = Xorshifts.Xoroshiro128Plus(1234)

@testset "test for bootstrap error throwing" begin
    cube_sum(x,y) = x^3+y^3
    x = [1,2,3]
    y = [4,5]
    @test_throws ErrorException bootstrap(cube_sum, x, y)

    x = [1,2,3]
    y = [4,5,6]
    z = [7,8,9]
    @test_throws MethodError bootstrap(cube_sum, x, y, z)

    cube_each(x,y) = [x^3, y^3]
    @test_throws ErrorException bootstrap(cube_each, x, y)
end

@testset "functional bootstrap test" begin
    mean_val =  5.0
    x =  fill(mean_val, 100) + randn!(rng, zeros(100))
    bs_error = Measurements.uncertainty(bootstrap(x))
    #check bootstraped std < original sample std
    @test bs_error < std(x)

    add(x,y,z) = x+y+z
    num_measurements_list = [10,100,1000,10000,100000]
    #create a set of normally random values about a mean value for x,y,z for different measuremets
    xs_list = [[fill(mean_val, num_measurements_list[n]) + randn!(rng, zeros(num_measurements_list[n])) for v in 1:3] for n in 1:length(num_measurements_list)]
    bootstrap_errors = [Measurements.uncertainty(bootstrap(add, xs_list[n][1], xs_list[n][2], xs_list[n][3])) for n in 1:length(num_measurements_list)]
    #test that error decreases as the number of observations gets larger
    @test issorted(bootstrap_errors, lt=Base.isgreater)
end

