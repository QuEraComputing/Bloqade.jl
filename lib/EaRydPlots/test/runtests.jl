using Test
using CairoMakie
using EaRydCore
using EaRydPlots

@testset "test plot" begin
    r = rand_state(12)
    plt = bitstring_histgram(r)
    @test size(plt.layout) == (1, 1)
end

r = rand_state(12)
plt = bitstring_histgram(r)
