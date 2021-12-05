using Test
using EaRydWaveforms
using EaRydWaveforms: edge_pad

wf = piecewise_linear(clocks=[0.0, 2.0, 3.0, 4.0], values=[0.0, 3.0, 1.1, 2.2])

using LinearAlgebra
clocks, values = [0.0, 1.0, 4.0, 5.0], [0.0, 3.0, 3.0, 0.0]
pad_clocks, pad_values = edge_pad(clocks, 8), edge_pad(values, 8)
f(x, b) = smooth(x, clocks, values, b) do t
    exp(-t^2/2)
end

function moving_f(x, b, n::Int)
    
    smooth(x, clocks, values)
end

xs = 0:1e-2:5
ys = [f(x, 0.3) for x in xs]

using GR
plot(xs, ys, backgroundcolor=255)

f(5.0, 0.1)
f(0.0, 0.1)



edge_pad([1, 2, 3, 4, 5], 2)

xs = rand(5)
ys = rand(5)
DSP.conv(xs, ys)

xs .* ys

simple_conv(xs, ys) = xs .* reverse(ys)
simple_conv(xs, ys)

a = [1, 2, 1, 2]
b = [1, 2, 3]
expectation = [1, 4, 8, 10, 7, 6]
DSP.conv(a, b)

xs = [0, 0, 1, 2, 1, 2, 0, 0]
ys = [3, 2, 1]
[sum(xs[idx:idx+2] .* ys) for idx in 1:6]
