using BloqadeSchema
using Bloqade
using JSON

# Maximal Rabi frequency and maximal (and minimal) detuning to see the effect of driving off-resonantly

atoms = AtomList([(0.0, 0.0)])

T=2
Ω = constant(;duration=T, value=8π)

# maximal detuning 
Δ = constant(;duration=T, value=40π)

# minimal detuning 
# Δ = constant(;duration=T, value=-40π)

ϕ = constant(;duration=T, value=0)

h = rydberg_h(atoms;Ω=Ω, Δ=Δ, ϕ=ϕ)
h_hardware, info = hardware_transform(h)
h = to_json(h_hardware)


open("lib/BloqadeSchema/schema_examples/aws_science_tasks/global_rydberg/global3.json","w") do f
    JSON.print(f, JSON.parse(h))
end


