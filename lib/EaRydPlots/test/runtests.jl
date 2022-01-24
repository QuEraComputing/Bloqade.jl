using Test
using CairoMakie
using EaRydCore
using EaRydPlots

@testset "test plot" begin
    r = rand_state(12)
    plt = bitstring_histgram(r)
    @test size(plt.layout) == (1, 1)
end

# r = rand_state(12)
# plt = bitstring_histgram(r)


# using EaRydWaveforms
# wf = piecewise_linear(clocks=[0.0, 1.0, 2.0, 3.0], values=[0.0, 2.0, 2.0, 1.0])

# EaRydPlots.draw(wf)
