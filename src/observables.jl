"""
    rydberg_density(reg, i::Int) -> Real

Calculates the rydberg density at site `i`.

```math
\\langle Op.n_i \\rangle
```
"""
rydberg_density(reg, i::Int) = real(expect(put(nqubits(reg), i=>Op.n), reg))

"""
    rydberg_density(reg) -> Vector

Return the rydberg density at each site.
"""
rydberg_density(reg) = [rydberg_density(reg, i) for i in 1:nqubits(reg)]

"""
    rydberg_corr([op=Op.n], reg) -> Matrix

Calculates the rydberg correlation.

```math
\\sum_{ij} \\langle op_i op_j \\rangle
```

# Arguments

- `op`: the correlation function, default is `Op.n`.
- `reg`: required, the register object.
"""
rydberg_corr(reg) = rydberg_corr(Op.n, reg)

function rydberg_corr(op, reg)
    return [
        expect(chain(nqubits(reg), put(i=>op), put(j=>op)), reg)
        for i in 1:nqubits(reg), j in 1:nqubits(reg)
    ]
end
