using Test
using Bloqade

reg = rand_state(10)
@test rydberg_density(reg, 2) ≈ expect(put(10, 2=>Op.n), reg)
@test rydberg_density(reg) ≈ [rydberg_density(reg, i) for i in 1:10]

@test rydberg_corr(reg) ≈ [
    expect(chain(nqubits(reg), put(i=>Op.n), put(j=>Op.n)), reg)
        for i in 1:nqubits(reg), j in 1:nqubits(reg)
]
