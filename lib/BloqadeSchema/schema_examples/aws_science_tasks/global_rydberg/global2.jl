using BloqadeSchema
using Bloqade
using JSON


# No detuning, and using a rabi_frequency pulse strength for a longer time, and see how much coherence is lost due to longer evolution time.

atoms = AtomList([(0.0, 0.0)])

T=2
epsilon =0.01
Ω = piecewise_constant(clocks=[0.0, epsilon, T+ epsilon , T+2*epsilon], values= 2π*[0.0, 3, 0])

H = rydberg_h(atoms;Ω=Ω)
h = to_json(H,waveform_tolerance=1e-1,warn=true)


open("lib/BloqadeSchema/schema_examples/aws_science_tasks/global_rydberg/global2.json","w") do f
    JSON.print(f, h)
end