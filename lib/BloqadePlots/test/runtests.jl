using Test
using BloqadePlots
using BloqadeWaveforms
using YaoArrayRegister
using YaoSubspaceArrayReg
using Random

wf = piecewise_linear(clocks=[0.0, 1.0, 2.0, 3.0], values=[0.0, 2.0, 2.0, 1.0])
draw(wf)

r = rand_state(12)
bitstring_hist(r; nlargest=10)

space = Subspace(5, randperm(1<<5)[1:10].-1)
r = rand_state(space)
bitstring_hist(r; nlargest=10)
