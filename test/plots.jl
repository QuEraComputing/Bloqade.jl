using Test
using Bloqade
using Random

wf = piecewise_linear(clocks = [0.0, 1.0, 2.0, 3.0], values = [0.0, 2.0, 2.0, 1.0])
Bloqade.plot(wf)

r = rand_state(12)
bitstring_hist(r; nlargest = 10)

space = Subspace(5, randperm(1 << 5)[1:10] .- 1)
r = rand_state(space)
bitstring_hist(r; nlargest = 10)

@testset "plot_densities" begin
    # dimension mismatch test
    atoms = generate_sites(SquareLattice(), 4, 4, scale=6.7)
    # Densities from adiabatic sweep on 3x3 Square Lattice
    densities = [0.9721147069025007,
                 0.010593663461943983,
                 0.9721147069025009,
                 0.010593663461943983,
                 0.987148044531143,
                 0.010593663461943978,
                 0.9721147069025005,
                 0.010593663461943981,
                 0.9721147069025007]
    @test_throws ArgumentError plot_densities(atoms, densities)
end
