"""
    rydberg_density(reg, i::Int)

Calculates the rydberg density at site `i`.

```math
\\langle Op.n_i \\rangle
```
"""
rydberg_density(reg, i::Int) = expect(put(nqubits(reg), i=>Op.n), reg)

"""
    rydberg_corr([op=Op.n], reg)

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
        expect(chain(put(i=>op), put(j=>op)), reg)
        for i in 1:nqubits(reg), j in 1:nqubits(reg)
    ]
end
