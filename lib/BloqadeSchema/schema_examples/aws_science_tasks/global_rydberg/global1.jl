using BloqadeSchema
using Bloqade
using JSON


# No detuning, and constant rabi_frequency pulse applied at max strength for a short period of time (to achieve a quick pi-pulse)


atoms = AtomList([(0.0, 0.0)])

T=1/8
epsilon =0.01
Ω = piecewise_constant(clocks=[0.0, epsilon, T+ epsilon , T+2*epsilon], values= 2π*[0.0, 4, 0])

H = rydberg_h(atoms;Ω=Ω)
h = to_json(H,waveform_tolerance=1e-1,warn=true)


open("lib/BloqadeSchema/schema_examples/aws_science_tasks/global_rydberg/global1.json","w") do f
    JSON.print(f, h)
end

