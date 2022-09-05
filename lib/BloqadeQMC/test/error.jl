using Test, BloqadeQMC
using Statistics: mean, std
using Measurements: value, uncertainty
using BloqadeQMC: jackknife
using Jackknife: estimate, variance, bias

@testset "simple numerical tests for jackknife" begin
    x = [3,2,11,7,6]
    
    # Resampling should not change estimate when f = mean. 
    m = mean(x)
    @test value(jackknife(mean, x)) ≈ m  

    # Resampling should give a better estimate for e.g. f = inv. We check this using the Julia library Jackknife
    g(x) = 1/mean(x) # The estimate function from Jackknife takes needs g to be explicity defined as a function of the mean
    @test value(jackknife(inv, x)) ≈ estimate(g, x)
    @test uncertainty(jackknife(inv, x)) ≈ sqrt(variance(g, x))

    # Test jackknife when dispatched on f = identity
    @test std(x) > uncertainty(jackknife(x))
end
