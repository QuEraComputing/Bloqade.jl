using Test
using EaRydWaveforms

clocks = collect(0:0.1:1)
values = cos.(clocks)
waveform = InterpolatedWaveform(clocks, values)
