"""
    unit_disk_graph(atoms::AbstractVector{<:RydAtom}, radius=1)

Create a unit disk graph from atom positions `atoms`. It returns a `Graphs.SimpleGraph` instance.
* `atoms` is vector of atoms positions.
* `radius` is the unit in the unit disk graph definition.
"""
function unit_disk_graph(atoms::AbstractVector{<:RydAtom}, radius=1)
    n_atoms = length(atoms)
    g = SimpleGraph(n_atoms)

    for k in 1:n_atoms
        indices = findall(atoms) do x
            distance(atoms[k], x) <= radius
        end

        for p in indices
            if k != p
                add_edge!(g, k, p)
            end
        end
    end

    return g
end

"""
    rand_unit_disk_graph(natoms::Int, ρ::Real)

Generate a random disk graph.

* `ρ` is defined as `n`/L^2, where L is the box size.
"""
function rand_unit_disk_graph(n::Int, ρ::Real)
    atoms = rand_atoms(n, ρ; ndims=2)
    unit_disk_graph(atoms, 1.0)
end
