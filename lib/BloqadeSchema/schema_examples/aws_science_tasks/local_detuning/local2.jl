using BloqadeSchema
using Bloqade
using JSON



# 10-50 atoms, separated by large distances (10um) to be outside each other blockade radius are excited with a pi-pulse. Gloal detuning is set to be negative, and local detuning is set so some atoms are on resonant with the pi-pulse, where others have resonance frequencies below or above the driving frequency.

a = 10;
N = 10;
atoms = generate_sites(ChainLattice(), N, scale=a)

# pi pulse for \Omega 
T=1/8
epsilon =0.01
Ωt = piecewise_linear(clocks=[0.0, epsilon, T+ epsilon , T+2*epsilon], values= 2π*[0.0, 4, 4,  0])

Δ_global = piecewise_linear(clocks=[0.0, epsilon, T+ epsilon , T+2*epsilon], values= 2π*[0.0,  -9, -9,  0])
Δ_local = piecewise_linear(clocks=[0.0, epsilon, T+ epsilon , T+2*epsilon], values= 2π*[0.0, 9, 9,  0])


Δt = map(1:length(atoms)) do idx
    if idx == floor(Int, N/2)+1
        Δ_global + 1* Δ_local
    else
        Δ_global + 0.1* Δ_local
    end
end ;

Δt = Vector{Waveform}(Δt)
H = rydberg_h(atoms; Δ = Δt, Ω = Ωt)
h = to_json(H, waveform_tolerance=1e-1, warn=true)

open("lib/BloqadeSchema/schema_examples/aws_science_tasks/local_detuning/local2.json","w") do f
    JSON.print(f, JSON.parse(h))
end







