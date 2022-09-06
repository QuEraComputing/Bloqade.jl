using BloqadeSchema
using Bloqade
using JSON



# preparing states with 1 detuning pattern. For more details, see https://queracomputing.github.io/Bloqade.jl/dev/tutorials/4.LGT/main/

a = 5.5;
N = 21;
atoms = generate_sites(ChainLattice(), N, scale=a);

total_time = 3.5;
Ωmax = 2π * 5;
Δmin = -2π * 10;
Δmax = 2π * 10;

Δq = 0.0;
tq = π/Ωmax;

Δ1 = piecewise_linear(clocks = [0.0, 0.2, total_time], values = [Δmin, Δmin, Δmax]);
Ω1 = piecewise_linear(clocks = [0.0, 0.2, total_time], values = [0.0, Ωmax, Ωmax]);


Δ2_single_defect = map(1:length(atoms)) do idx
    if idx == floor(Int, N/2)+1
        append(Δ1, constant(duration=tq, value=Δq))
    else
        append(Δ1, constant(duration=tq, value=Δmax))
    end
end ;

Ω2 = append(Ω1, constant(duration=tq, value=Ωmax));

H = rydberg_h(atoms; Δ = Δ2_single_defect, Ω = Ω2)
h = BloqadeSchema.to_json(H,waveform_tolerance=1e-1, warn=true)

open("generate_schema/1detuning_waveform.json","w") do f
    JSON.print(f, h)
end

# test if the local detuning is generated correctly 

h2 = BloqadeSchema.from_json(h)
(atoms,phi,Omega,Delta) = BloqadeSchema.get_rydberg_params(h2)
