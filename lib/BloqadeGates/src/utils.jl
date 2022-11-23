"""
    two_level_indices(n::Integer)

Return the indices of the `n` qubits in the `n`-site 3-level system. 
It is used to test the correctness of the pulse sequence generated for the given gate.
"""
function two_level_indices(n::Integer)
    n == 1 && return [1, 2]
    return [two_level_indices(n-1); two_level_indices(n-1).+3^(n-1)]
end