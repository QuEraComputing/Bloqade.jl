using BloqadeSchema
using Bloqade
using JSON



# preparing 1D Z2 phase. For more details, see https://queracomputing.github.io/Bloqade.jl/dev/tutorials/2.adiabatic/main/

total_time = 3.0;
Ω_max = 2π * 4;
Ω = piecewise_linear(clocks = [0.0, 0.1, 2.1, 2.2, total_time], values = [0.0, Ω_max, Ω_max, 0, 0]);

U1 = -2π * 10;
U2 = 2π * 10;
Δ = piecewise_linear(clocks = [0.0, 0.6, 2.1, total_time], values = [U1, U1, U2, U2]);

ϕ = constant(;duration=total_time, value=0)

nsites = 9
atoms = generate_sites(ChainLattice(), nsites, scale = 5.72)


h = rydberg_h(atoms; Δ, Ω, ϕ)
h_hardware, info = hardware_transform(h)
h = to_json(h_hardware)


open("lib/BloqadeSchema/schema_examples/tutorial_examples/1D_Z2.json","w") do f
    JSON.print(f, JSON.parse(h))
end