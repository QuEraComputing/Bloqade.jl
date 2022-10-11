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

ϕ = constant(duration = 4.2, value = 0);

h = rydberg_h(atoms; Δ = Δ_tot, Ω = Ω_tot, ϕ=ϕ)
device_capabilities = get_device_capabilities()
# WARNING: This will result in invalid JSON file for actual QPU.
device_capabilities.rydberg.global_value.timeMax = 5
h_hardware, info = hardware_transform(h, device_capabilities=device_capabilities)
h = to_json(h_hardware; device_capabilities)

open("lib/BloqadeSchema/schema_examples/tutorial_examples/1D_scar.json","w") do f
    JSON.print(f, JSON.parse(h))
end



