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

# Throw an error if there are no waveforms supplied for Ω, Δ, or ϕ
@test_throws ErrorException get_average_rydberg_densities(
                              generate_sites(ChainLattice(), 1), 
                              zero_state(1);
                              dt=0.1)

# Throw an error if waveforms are present but durations are incompatible
## Waveforms taken from "Adiabatic Evolution" Tutorial of Bloqade Documentation: 
## https://queracomputing.github.io/Bloqade.jl/dev/tutorials/2.adiabatic/main/#Preparation-of-Ordered-States-in-1D

Ω_max = 2π * 4;
Ω = piecewise_linear(clocks = [0.0, 0.1, 1.2, 2.2, 3.0], values = [0.0, Ω_max, Ω_max, 0, 0 ]) 

U1 = -2π * 10;
U2 = 2π * 10;
Δ = piecewise_linear(clocks = [0.0, 0.6, 2.1, 4.1],values = [U1, U1, U2, U2]);
@test_throws ErrorException get_average_rydberg_densities(
                              generate_sites(ChainLattice(), 1), 
                              zero_state(1);
                              Ω = Ω,
                              Δ = Δ,
                              dt=0.1)