using BloqadeSchema
using Bloqade
using JSON

# Maximal Rabi frequency and maximal (and minimal) detuning to see the effect of driving off-resonantly

atoms = AtomList([(0.0, 0.0)])

T=2
epsilon =0.01
Ω = piecewise_constant(clocks=[0.0, epsilon, T+ epsilon , T+2*epsilon], values= 2π*[0.0, 4, 0])

# maximal detuning 
Δ = piecewise_constant(clocks=[0.0, epsilon, T+ epsilon , T+2*epsilon], values= 2π*[0.0, 20, 0])

# minimal detuning 
# Δ = piecewise_constant(clocks=[0.0, epsilon, T+ epsilon , T+2*epsilon], values= 2π*[0.0, -20, 0])

H = rydberg_h(atoms;Ω=Ω, Δ = Δ)
h = BloqadeSchema.to_json(H,waveform_tolerance=1e-1,warn=true)


open("generate_schema/aws_science_tasks/global_rydberg/global3.json","w") do f
    JSON.print(f, h)
end



