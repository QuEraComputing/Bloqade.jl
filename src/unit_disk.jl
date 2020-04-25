using StatsBase
export AtomPosition, unit_disk_graph

struct AtomPosition
    coordinates::Matrix{Float64}
end

"""
    AtomPosition(n::Int, ndims::Int=2)

Create a random atom position of `n` atoms in `ndims` space.
"""
AtomPosition(n::Int, ndims::Int=2) = AtomPosition(randn(ndims, n))

"""
    AtomPosition(n::Int, ρ::Float64, ndims::Int=2)

Create a random atom position of `n` atoms and with density `ρ` in `ndims` space.
The size of the box is L^ndims, where the linear dimension L = (n/ρ)^(1/ndims).
"""
function AtomPosition(n::Int, ρ::Float64, ndims::Int=2)
    L = (n/ρ)^(1/ndims)
    atom_coordinates = rand(ndims, n) * L
    AtomPosition(atom_coordinates)
end

"""
    AtomPosition(n::Int, ff::Float64, geometry::String)

Create a random atom position of `n` atoms and with a filling factor `ff`
in 2D space on a `square` lattice.
"""
function AtomPosition(n::Int, ff::Float64, geometry::String)
    if geometry == "square"
        L = round(Int64,sqrt(n/ff))
        atom_coordinates_linear = sample(1:L^2,n,replace = false)
        atom_coordinates_x = (atom_coordinates_linear .- 1) .÷ L .+ 1
        atom_coordinates_y = (atom_coordinates_linear .- 1) .% L .+ 1
        atom_coordinates = vcat(atom_coordinates_x', atom_coordinates_y')
        AtomPosition(atom_coordinates)
    end
end

Base.ndims(x::AtomPosition) = size(x.coordinates, 1)
Base.length(x::AtomPosition) = size(x.coordinates, 2)
Base.getindex(x::AtomPosition, k::Int) = getindex(x.coordinates, :, k)
Base.keys(x::AtomPosition) = LinearIndices(1:length(x))

function Base.iterate(it::AtomPosition, st=1)
    if st > length(it)
        return
    else
        return it[st], st + 1
    end
end

function unit_disk_graph(atoms::AtomPosition, radius::Float64=1)
    n_atoms = length(atoms)
    g = SimpleGraph(n_atoms)

    for k in 1:n_atoms
        indices = findall(atoms) do x
            norm(atoms[k] - x) <= radius
        end

        for p in indices
            if k != p
                add_edge!(g, k, p)
            end
        end
    end
    return g
end
