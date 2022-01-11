using Makie
using CairoMakie
using MakieCore
using EaRydCore

@macroexpand @recipe(BitStringHistgram) do scene
end


r = rand_state(12)
plt = bitstring_histgram(r)
