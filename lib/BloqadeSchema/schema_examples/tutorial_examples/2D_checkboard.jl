using BloqadeSchema
using Bloqade
using JSON


# preparing 2D checkboard phase. For more details, see https://queracomputing.github.io/Bloqade.jl/dev/tutorials/2.adiabatic/main/

nx, ny = 3, 3
nsites = nx * ny
atoms = generate_sites(SquareLattice(), nx, ny, scale = 6.7)

total_time = 2.9
Ω_max = 2π * 4.3
Ω = piecewise_linear(clocks = [0.0, 0.3, 2.6, total_time], values = [0.0, Ω_max, Ω_max, 0]);

U = 2π * 15.0
Δ = piecewise_linear(clocks = [0.0, 0.3, 2.6, total_time], values = [-U, -U, U, U]);

H = rydberg_h(atoms; Δ, Ω)
h = to_json(H,waveform_tolerance=1e-1, warn=true)

open("lib/BloqadeSchema/schema_examples/tutorial_examples/2D_checkboard.json","w") do f
    JSON.print(f, h)
end
