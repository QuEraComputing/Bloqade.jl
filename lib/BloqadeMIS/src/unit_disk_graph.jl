"""
    unit_disk_graph(atoms::AbstractVector, radius=1)

Create a unit disk graph from atom positions `atoms`. It returns a `Graphs.SimpleGraph` instance.
* `atoms` is vector of atoms positions.
* `radius` is the unit in the unit disk graph definition.
"""
function unit_disk_graph(atoms::AbstractVector, radius = 1)
    n_atoms = length(atoms)
    g = SimpleGraph(n_atoms)

    for k in 1:n_atoms
        indices = findall(atoms) do x
            return BloqadeExpr.distance(atoms[k], x) <= radius
        end

        for p in indices
            if k != p
                add_edge!(g, k, p)
            end
        end
    end

    return g
end
