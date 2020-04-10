export AtomPosition, unit_disk_graph

struct AtomPosition
    coordinates::Matrix{Float64}
end

"""
    AtomPosition(n::Int, ndims::Int=2)

Create a random atom position of `n` atoms in `ndims` space.
"""
AtomPosition(n::Int, ndims::Int=2) = AtomPosition(randn(ndims, n))

Base.ndims(x::AtomPosition) = size(x.coordinates, 1)
Base.length(x::AtomPosition) = size(x.coordinates, 2)
Base.getindex(x::AtomPosition, k::Int) = getindex(x.coordinates, :, k)
Base.keys(x::AtomPosition) = LinearIndices(1:length(x))

function Base.iterate(it::AtomPosition, st=1)
    if st > size(it.coordinates, 2)
        return
    else
        return it.coordinates[:, st], st + 1
    end
end

function unit_disk_graph(atoms::AtomPosition, radius::Float64)
    n_atoms = size(atoms.coordinates, 2)
    g = SimpleGraph(n_atoms)

    for k in 1:n_atoms
        indices = findall(atoms) do x
            norm(atoms.coordinates[:, k] - x) < radius
        end

        for p in indices
            if k != p
                add_edge!(g, k, p)
            end
        end
    end
    return g
end
