using Test, BloqadeQMC
using Statistics: mean, std
using Measurements: value, uncertainty
using BloqadeQMC: jackknife
using Random
using RandomNumbers

rng = Xorshifts.Xoroshiro128Plus(1234)

@testset "simple numerical tests for jackknife" begin
    x = [3,2,11,7,6]
    
    # Resampling should not change estimate when f = mean. 
    m = mean(x)
    @test value(jackknife(mean, x)) ≈ m  

    # Resampling should give a better estimate for e.g. f = inv. We check this using the Julia library Jackknife
    @test value(jackknife(inv, x)) ≈ 0.158110765
    @test uncertainty(jackknife(inv, x)) ≈ 0.052468798
end

@testset "functional tests for jackknife" begin

    mean_val =  5.0
    x =  fill(mean_val, 100) + randn!(rng, zeros(100))
    j_error = uncertainty(jackknife(x))
    #check jackknifed std < original sample std
    @test j_error < std(x)  

    add(x,y,z) = x+y+z
    num_measurements_list = [10,100,1000,10000,100000]
    #create a set of normally random values about a mean value for x,y,z for different measuremensts
    xs_list = [[fill(mean_val, num_measurements_list[n]) + randn!(rng, zeros(num_measurements_list[n])) for v in 1:3] for n in 1:length(num_measurements_list)]
    jackknife_errors = [uncertainty(jackknife(add, xs_list[n][1], xs_list[n][2], xs_list[n][3])) for n in 1:length(num_measurements_list)]
    #test that error decreases as the number of observations gets larger
    @test issorted(jackknife_errors, lt=Base.isgreater)

end
