using BloqadeSchema
using Bloqade
using JSON



# Set both global and local detuning to their minimum and maximum values in different tasks to see if the machine can handle it and the dynamics is correct (see experiment 1))
# two atoms are far seperated with each other 

a = 11;
N = 2;
atoms = generate_sites(ChainLattice(), N, scale=a);

total_time = 3.5;
Ωmax = 2π * 4;
Δmin = -2π * 10;
Δmax = 2π * 10;


Δ_global = piecewise_linear(clocks = [0.0, 0.2, total_time], values = [Δmin, Δmin, Δmax])
Δ_local = piecewise_linear(clocks = [0.0, 0.9, total_time], values = [0, 1/2* Δmax , Δmax])

Δt = map(1:length(atoms)) do idx
    if idx == floor(Int, N/2)+1
        Δ_global + 1* Δ_local
    else
        Δ_global + 2* Δ_local
    end
end ;

Ωt = piecewise_linear(clocks = [0.0, 0.2, total_time], values = [0.0, Ωmax, Ωmax]);

H = rydberg_h(atoms; Δ = Δt , Ω = Ωt)
h = BloqadeSchema.to_json(H,waveform_tolerance=1e-1, warn=true)

open("generate_schema/aws_science_tasks/local_detuning/local1.json","w") do f
    JSON.print(f, h)
end
