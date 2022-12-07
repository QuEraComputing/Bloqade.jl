using Test
using BloqadeLattices
using BloqadeExpr
using BloqadeQMC: make_prob_vector, ImprovedOperatorSampler


@testset "rydberg_qmc" begin
    nsites = 3
    atoms = generate_sites(ChainLattice(), nsites, scale = 5.48)
    lattice = parallelepiped_region(ChainLattice(), (nsites,); scale = 5.48)

    Ω = 2π * 4
    Δ = 2π * 4
    Ω_N = Ω*ones(nsites)
    Δ_N = Δ*ones(nsites)

    for positions in [atoms, lattice]
        h = rydberg_h(positions; Δ, Ω)
        H = rydberg_qmc(h)
    
        V = rydberg_interaction_matrix(positions, h.rydberg_term.C)
    
        ops, p, energy_shift = make_prob_vector(AbstractRydberg, V, Ω_N, Δ_N, epsilon=0.0)
        op_sampler = ImprovedOperatorSampler(AbstractLTFIM, ops, p)
    
        @test H.op_sampler.operators == op_sampler.operators
        @test H.op_sampler.op_log_weights == op_sampler.op_log_weights
        @test H.V == V
        @test H.Ω == Ω_N
        @test H.δ == Δ_N
        @test H.atoms == positions
        @test H.energy_shift == energy_shift
    end

end

