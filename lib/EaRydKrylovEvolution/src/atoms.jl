"""
    AbstractAtom

Abstract type for atoms.
"""
abstract type AbstractAtom end

"""
    RydAtom{N,T} <: AbstractAtom

Position of one Rydberg atom.

    RydAtom(locs::NTuple{N, T})

Create a position of one Rydberg atom. The default unit is `μm`.
"""
struct RydAtom{N,T <: Real} <: AbstractAtom
    loc::NTuple{N,T}
end

RydAtom(locs::Tuple) = RydAtom(default_unit(μm, locs))

# constructors
"""
    Atom2D{T} = RydAtom{2,T}
"""
const Atom2D{T} = RydAtom{2,T}

"""
    RydAtom(locations...)

Create a `RydAtom` from given locations.
"""
RydAtom(args...) = RydAtom(args)

"""
    RydAtom(locations::Vector)

Create a `RydAtom` from given list of locations.
"""
RydAtom(x::AbstractVector) = RydAtom(x...)
# NOTE: the following code is not fast

"""
    rydatoms(::AbstractMatrix)

Create a list of [`RydAtom`](@ref)s from a nx2 location matrix.
"""
rydatoms(locs::AbstractMatrix) = [RydAtom(locs[:,i]...) for i=1:size(locs,2)]

# interfaces
Base.ndims(x::RydAtom{N}) where N = N
Base.getindex(x::RydAtom, k::Int) = getindex(x.loc, k)
Base.length(::RydAtom{N}) where N = N
Base.eltype(::RydAtom{N, T}) where {N, T} = T
Base.iterate(x::RydAtom, st...) = iterate(x.loc, st...)
axis(a::RydAtom) = a.loc

"""
    distance(a::RydAtom, b::RydAtom)

Return the distance between two Rydberg atoms.
"""
distance(a::RydAtom, b::RydAtom) = sqrt(mapreduce(x->x^2, +, axis(a) .- axis(b)))

"""
    rand_atoms(n::Int, ρ::Float64; ndims::Int=2)

Create a random atom position of `n` atoms and with density `ρ` in `ndims` space.
The size of the box is ``L^ndims``, where the linear dimension ``L = (n/ρ)^(1/ndims)``.
"""
function rand_atoms(n::Int, ρ::Float64; ndims::Int=2)
    L = (n/ρ)^(1/ndims)
    atoms = rydatoms(rand(ndims, n) .* L)
    sort_atoms!(atoms)
    return SVector{length(atoms)}(atoms)
end

"""
    square_lattice(n::Int, ff::Float64)

Create a list of Rydberg atoms on a square lattice with given size `n` and
filling factor `ff`.
"""
function square_lattice(n::Int, ff::Float64)
    L = ceil(Int64,sqrt(n/ff))
    atom_coordinates_linear = sample(1:L^2,n,replace = false)
    atom_coordinates_x = (atom_coordinates_linear .- 1) .÷ L .+ 1
    atom_coordinates_y = (atom_coordinates_linear .- 1) .% L .+ 1
    atom_coordinates = vcat(atom_coordinates_x', atom_coordinates_y')
    atoms = sort_atoms!(rydatoms(atom_coordinates))
    return SVector{length(atoms)}(atoms)
end

"""
    sort_atoms!(atoms; by=nothing)

Sort a list of atom positions. Sort by coordinates if `by` is `nothing`.
"""
function sort_atoms!(atoms::Vector{<:RydAtom}; by=nothing)
    if isnothing(by)
        sort!(atoms, by=p->p.loc)
    else
        sort!(atoms, by=by)
    end
    return atoms
end
