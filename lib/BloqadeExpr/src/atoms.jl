# convenient mapping functions for atoms

"""
    matrix_to_positions(locs::AbstractMatrix{T}) -> Vector{NTuple{D, T}}

Convert a dxn location matrix to a list of positions.
"""
function matrix_to_positions(locs::AbstractMatrix)
    return [(locs[:,i]..., ) for i in 1:size(locs, 2)]
end
