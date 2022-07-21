using Test
using Bloqade

reg = rand_state(10)
@test rydberg_density(reg, 2) ≈ expect(put(10, 2 => Op.n), reg)
@test rydberg_density(reg) ≈ [rydberg_density(reg, i) for i in 1:10]

@test rydberg_corr(reg) ≈
      [expect(chain(nqubits(reg), put(i => Op.n), put(j => Op.n)), reg) for i in 1:nqubits(reg), j in 1:nqubits(reg)]


# Test 1 for get_average_rydberg_densities: π-pulse for 1-atom
@test [[0.0], [1.0]] ≈ get_average_rydberg_densities(
      generate_sites(ChainLattice(), 1), 
      zero_state(1); 
      Ω = piecewise_linear(clocks = [0.0, 0.1], values = [π/(0.1), π/(0.1)]), 
      dt=0.1)
