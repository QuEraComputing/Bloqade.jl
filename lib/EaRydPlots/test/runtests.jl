using Test
using EaRydPlots
using EaRydWaveforms
using YaoArrayRegister
using YaoSubspaceArrayReg
using Random

wf = piecewise_linear(clocks=[0.0, 1.0, 2.0, 3.0], values=[0.0, 2.0, 2.0, 1.0])
fig = plt.figure()
draw(wf)
fig


r = rand_state(12)
fig = plt.figure()
bitstring_hist(r; nlargest=10)
fig

space = Subspace(5, randperm(1<<5)[1:10].-1)
r = rand_state(space)
fig = plt.figure()
bitstring_hist(r; nlargest=10)
fig
