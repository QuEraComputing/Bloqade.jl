using BloqadeSchema
using Bloqade
using JSON


#  1D scar dynamics. For more details, see https://queracomputing.github.io/Bloqade.jl/dev/tutorials/3.quantum-scar/main/

nsites = 9
atoms = generate_sites(ChainLattice(), nsites, scale = 5.72)

Δ1 = piecewise_linear(clocks = [0.0, 0.3, 1.6, 2.2], values = 2π * [-10.0, -10.0, 10.0, 10.0]);
Ω1 = piecewise_linear(clocks = [0.0, 0.05, 1.6, 2.2], values = 2π * [0.0, 4.0, 4.0, 0.0]);

Ω2 = constant(duration = 2.0, value = 2 * 2π);
Δ2 = constant(duration = 2.0, value = 0);

Ω_tot = append(Ω1, Ω2);
Δ_tot = append(Δ1, Δ2);

H = rydberg_h(atoms; Δ = Δ_tot, Ω = Ω_tot)
h = BloqadeSchema.to_json(H,waveform_tolerance=1e-1, warn=true)

open("generate_schema/tutorial_examples/1D_scar.json","w") do f
    JSON.print(f, h)
end



