using EaRyd
using CairoMakie
using EaRydPlots

nx, ny = 3, 3
nsites = nx*ny
# pulse 
Ω_max = 2.3 * 2 * pi
Ω = piecewise_linear(clocks=[0.0, 0.252, 1.052, 1.6], values=[0.0, Ω_max , Ω_max , 0])

U = Ω_max / 2.3
Δ = piecewise_linear(clocks=[0.0, 0.252, 1.052, 1.6], values=[-6*U, -6*U, 2*U , 2*U])

# lattice
atoms = generate_sites(SquareLattice(), nx, ny, scale=9.629)

h = rydberg_h(atoms; C=2 * pi * 858386, Δ, Ω)
prob = ODEEvolution(zero_state(9), 1.6, h)
emulate!(prob) # run the time evolution directly

densities = map(1:nsites) do i
    real(expect(put(nsites, i=>Op.n), prob.reg))
end

bitstring_histgram(prob.reg; nlargest_possible=20)
