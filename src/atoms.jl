abstract type AbstractAtom end

"""
    RydAtom{N,T} <: AbstractAtom
    RydAtom(locations...)
    RydAtom(locations)

Rydberg atom.
"""
struct RydAtom{N,T} <: AbstractAtom
    loc::NTuple{N,T}
end

# constructors
"""
    Atom2D{T} = RydAtom{2,T}
"""
const Atom2D{T} = RydAtom{2,T}
RydAtom(args...) = RydAtom(args)
RydAtom(x::AbstractVector) = RydAtom(x...)
# NOTE: the following code is not fast
rydatoms(locs::AbstractMatrix) = [RydAtom(locs[:,i]...) for i=1:size(locs,2)]

# interfaces
Base.ndims(x::RydAtom{N}) where N = N
Base.getindex(x::RydAtom, k::Int) = getindex(x.loc, k)

axis(a::RydAtom) = a.loc
distance(a::RydAtom, b::RydAtom) = sqrt(mapreduce(x->x^2, +, axis(a) .- axis(b)))

"""
    rand_atoms(n::Int, ρ::Float64; ndims::Int=2)

Create a random atom position of `n` atoms and with density `ρ` in `ndims` space.
The size of the box is ``L^ndims``, where the linear dimension ``L = (n/ρ)^(1/ndims)``.
"""
function rand_atoms(n::Int, ρ::Float64; ndims::Int=2)
    L = (n/ρ)^(1/ndims)
    rydatoms(rand(ndims, n) .* L)
end

"""
    lattice_atoms(n::Int, ff::Float64, geometry::String)

Create a random atom position of `n` atoms and with a filling factor `ff`
in 2D space on a `square` lattice.
"""
function lattice_atoms(n::Int, ff::Float64, geometry=:square)
    if geometry == :square
        square_lattice(n, ff)
    end
end

function square_lattice(n::Int, ff::Float64)
    L = ceil(Int64,sqrt(n/ff))
    atom_coordinates_linear = sample(1:L^2,n,replace = false)
    atom_coordinates_x = (atom_coordinates_linear .- 1) .÷ L .+ 1
    atom_coordinates_y = (atom_coordinates_linear .- 1) .% L .+ 1
    atom_coordinates = vcat(atom_coordinates_x', atom_coordinates_y')
    rydatoms(atom_coordinates)
end
