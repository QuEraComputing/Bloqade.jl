"""
    AbstractLattice{D}

Supertype for all `D` dimensional lattices.

# Implementation

[`lattice_vectors`](@ref) and [`lattice_sites`](@ref) functions must be defined
which should both return an indexable iterable containing the Bravais lattice vectors and 
lattice sites respectively. (e.g.: [`GeneralLattice`](@ref) returns a tuple of tuples containing the 
Bravais lattice vectors and lattice sites).
"""
abstract type AbstractLattice{D} end

"""
    dimension(::AbstractLattice{D})

Returns the space dimension of target lattice.
e.g. [`ChainLattice`](@ref) is a 1D lattice, hence returns 1.
"""
dimension(::AbstractLattice{D}) where {D} = D

function _generate_sites(lattice_vectors, lattice_sites, repeats::Vararg{Int,D}; scale = 1.0) where {D}
    @assert length(lattice_vectors) == D
    @assert D > 0 && length(lattice_sites) > 0
    @assert all(>=(0), repeats)
    @assert all(x -> length(x) == (D), lattice_vectors)
    @assert all(x -> length(x) == (D), lattice_sites)
    T = eltype(first(lattice_vectors))
    locations = NTuple{D,T}[]  # we might want to avoid using `push!` later.
    for ci in CartesianIndices(repeats)
        baseloc = mapreduce(i -> (ci.I[i] - 1) .* lattice_vectors[i], (x, y) -> x .+ y, 1:D)
        for siteloc in lattice_sites
            push!(locations, (baseloc .+ siteloc) .* scale)
        end
    end
    return locations
end

"""
    GeneralLattice{D,K,T} <: AbstractLattice{D}
    GeneralLattice(vectors, sites)

The general lattice type for tiling the space. Type parameter `D` is the dimension,
`K` is the number of sites in a unit cell and `T` is the data type for coordinates, e.g. `Float64`. Input arguments are

* `vectors` is a vector/tuple of D-tuple. Its length is D, it specifies the Bravais lattice vectors.
* `sites` is a vector/tuple of D-tuple. Its length is K, it specifies the sites inside a Bravais cell.
"""
struct GeneralLattice{D,K,T} <: AbstractLattice{D}
    vectors::NTuple{D,NTuple{D,T}}
    sites::NTuple{K,NTuple{D,T}}
end
GeneralLattice(vectors, sites) = GeneralLattice((vectors...,), (sites...,))

"""
    lattice_vectors(lattice::AbstractLattice)

Returns Bravais lattice vectors as a D-Tuple of D-Tuple, where D is the space dimension.
"""
lattice_vectors(general_lattice::GeneralLattice) = general_lattice.vectors

"""
    lattice_sites(lattice::AbstractLattice)

Returns sites in a Bravais lattice unit cell as a Tuple of D-Tuple, where D is the space dimension.
"""
lattice_sites(general_lattice::GeneralLattice) = general_lattice.sites

"""
    struct HoneycombLattice <: AbstractLattice{2}

Type representing 2D Honeycomb Lattice.

Overriden functions to return lattice vectors and sites exists as 
[`lattice_vectors(::HoneycombLattice)`](@ref) and
[`lattice_sites(::HoneycombLattice)`](@ref).

# Fields

None
"""
struct HoneycombLattice <: AbstractLattice{2} end

"""
    lattice_vectors(::HoneycombLattice)

Returns the Bravais lattice vectors for a Honeycomb lattice as a Tuple of Tuples containing
floats.

The vectors are defined as:
- 𝐚₁ = (1.0, 0.0)
- 𝐚₂ = (0.5, 0.5√3)
"""
lattice_vectors(::HoneycombLattice) = ((1.0, 0.0), (0.5, 0.5 * sqrt(3)))

"""
    lattice_sites(::HoneycombLattice)

Returns the Bravais Lattice sites for a Honeycomb lattice as a Tuple of Tuples containing
floats.

The sites are defined as:
- (0.0, 0.0)
- (0.5, 0.5√3)
"""
lattice_sites(::HoneycombLattice) = ((0.0, 0.0), (0.5, 0.5 / sqrt(3)))

"""
    struct SquareLattice <: AbstractLattice{2}

Type representing 2D Square Lattice.

Overriden functions to return lattice vectors and sites exists as 
[`lattice_vectors(::SquareLattice)`](@ref) and
[`lattice_sites(::SquareLattice)`](@ref).

# Fields

None
"""
struct SquareLattice <: AbstractLattice{2} end

"""
    lattice_vectors(::SquareLattice)

Returns the Bravais lattice vectors for a Square lattice as a Tuple of Tuples containing
floats.
    
The vectors are defined as:
- 𝐚₁ = (1.0, 0.0)
- 𝐚₂ = (0.0, 1.0)
"""
lattice_vectors(::SquareLattice) = ((1.0, 0.0), (0.0, 1.0))

"""
    lattice_sites(::SquareLattice)

Returns the Bravais Lattice sites for a Square lattice as a Tuple of Tuples containing
floats.

The sites are defined as:
- (0.0, 0.0)
"""
lattice_sites(::SquareLattice) = ((0.0, 0.0),)

"""
    struct TriangularLattice <: AbstractLattice{2}

Type representing 2D Square Lattice.

Overriden functions to return lattice vectors and sites exists as 
[`lattice_vectors(::TriangularLattice)`](@ref) and
[`lattice_sites(::TriangularLattice)`](@ref).

# Fields

None
"""
struct TriangularLattice <: AbstractLattice{2} end

"""
    lattice_vectors(::TriangularLattice)

Returns the Bravais lattice vectors for a Triangular lattice as a Tuple of Tuples containing
floats.
    
The vectors are defined as:
- 𝐚₁ = (1.0, 0.0)
- 𝐚₂ = (0.5, 0.5√3)
"""
lattice_vectors(::TriangularLattice) = ((1.0, 0.0), (0.5, 0.5 * sqrt(3)))

"""
    lattice_sites(::TriangularLattice)

Returns the Bravais Lattice sites for a Triangular lattice as a Tuple of Tuples containing
floats.

The sites are defined as:
- (0.0, 0.0)
"""
lattice_sites(::TriangularLattice) = ((0.0, 0.0),)

"""
    struct ChainLattice <: AbstractLattice{1}

Type representing 1D Chain Lattice.

Overriden functions to return lattice vectors and sites exists as 
[`lattice_vectors(::ChainLattice)`](@ref) and
[`lattice_sites(::ChainLattice)`](@ref).

# Fields

None
"""
struct ChainLattice <: AbstractLattice{1} end

"""
    lattice_vectors(::ChainLattice)

Returns the Bravais lattice vectors for a Chain lattice as a Tuple of Tuples containing
floats.
    
The vectors are defined as:
- 𝐚₁ = (1.0,)
"""
lattice_vectors(::ChainLattice) = ((1.0,),)

"""
    lattice_sites(::ChainLattice)

Returns the Bravais Lattice sites for a Chain lattice as a Tuple of Tuples containing
floats.

The sites are defined as:
- (0.0,)
"""
lattice_sites(::ChainLattice) = ((0.0,),)

"""
    struct LiebLattice <: AbstractLattice{2}

Type representing 1D Chain Lattice.

Overriden functions to return lattice vectors and sites exists as 
[`lattice_vectors(::LiebLattice)`](@ref) and
[`lattice_sites(::LiebLattice)`](@ref).

# Fields

None
"""
struct LiebLattice <: AbstractLattice{2} end

"""
    lattice_vectors(::LiebLattice)

Returns the Bravais lattice vectors for a Lieb lattice as a Tuple of Tuples containing
floats.
        
The vectors are defined as:
- 𝐚₁ = (1.0, 0.0)
- 𝐚₂ = (0.0, 1.0)
"""
lattice_vectors(::LiebLattice) = ((1.0, 0.0), (0.0, 1.0))

"""
    lattice_sites(::LiebLattice)

Returns the Bravais Lattice sites for a Lieb lattice as a Tuple of Tuples containing
floats.

The sites are defined as:
- (0.0, 0.0)
- (0.5, 0.0)
- (0.0, 0.5)
"""
lattice_sites(::LiebLattice) = ((0.0, 0.0), (0.5, 0.0), (0.0, 0.5))

"""
    struct KagomeLattice <: AbstractLattice{2}

Type representing 2D Kagome Lattice.

Overriden functions to return lattice vectors and sites exists as 
[`lattice_vectors(::KagomeLattice)`](@ref) and
[`lattice_sites(::KagomeLattice)`](@ref).

# Fields

None
"""
struct KagomeLattice <: AbstractLattice{2} end

"""
    lattice_vectors(::KagomeLattice)

Returns the Bravais lattice vectors for a Kagome lattice as a Tuple of Tuples containing
floats.
        
The vectors are defined as:
- 𝐚₁ = (1.0, 0.0)
- 𝐚₂ = (0.5, 0.5√3)
"""
lattice_vectors(::KagomeLattice) = ((1.0, 0.0), (0.5, 0.5 * sqrt(3)))

"""
    lattice_sites(::KagomeLattice)

Returns the Bravais Lattice sites for a Lieb lattice as a Tuple of Tuples containing
floats.

The sites are defined as:
- (0.0, 0.0)
- (0.25, 0.25√3)
- (0.75, 0.25√3)
"""
lattice_sites(::KagomeLattice) = ((0.0, 0.0), (0.25, 0.25 * sqrt(3)), (0.75, 0.25 * sqrt(3)))

"""
    struct RectangularLattice <: AbstractLattice{2}

Type representing 2D Rectangular Lattice.

Overriden functions to return lattice vectors and sites exists as 
[`lattice_vectors(::RectangularLattice)`](@ref) and
[`lattice_sites(::RectangularLattice)`](@ref).

# Fields
- `aspect_ratio::Float64`: defines the length of one of the Bravais lattice vectors
(see [`lattice_vectors(r::RectangularLattice)`](@ref))
"""
struct RectangularLattice <: AbstractLattice{2}
    aspect_ratio::Float64
end

"""
    lattice_vectors(r::RectangularLattice)

Returns the Bravais lattice vectors for a Rectangular lattice as a Tuple of Tuples containing
floats.
        
The vectors are defined as:
- 𝐚₁ = (1.0, 0.0)
- 𝐚₂ = (0.0, `r.aspect_ratio`), where `aspect_ratio` is a `Float64`.
"""
lattice_vectors(r::RectangularLattice) = ((1.0, 0.0), (0.0, r.aspect_ratio))

"""
    lattice_sites(::RectangularLattice)

Returns the Bravais Lattice sites for a Rectangular lattice as a Tuple of Tuples containing
floats.

The sites are defined as:
- (0.0, 0.0)
"""
lattice_sites(::RectangularLattice) = ((0.0, 0.0),)

"""
    AtomList{D, T} <: AbstractVector{NTuple{D, T}}
    AtomList(atoms::Vector{<:NTuple})

A list of atoms in `D` dimensional space.
"""
struct AtomList{D,T} <: AbstractVector{NTuple{D,T}}
    atoms::Vector{NTuple{D,T}}
end

Base.size(list::AtomList) = size(list.atoms)
Base.length(list::AtomList) = length(list.atoms)
Base.getindex(list::AtomList, idx::AbstractRange) = AtomList(list.atoms[idx])
Base.getindex(list::AtomList, idx::AbstractVector) = AtomList(list.atoms[idx])
Base.getindex(list::AtomList, idx::Int) = list.atoms[idx]

"""
    generate_sites(lattice::AbstractLattice{D}, repeats::Vararg{Int,D}; scale=1.0)

Returns an [`AtomList`](@ref) instance by tiling the specified `lattice`.
The tiling repeat the `sites` of the lattice `m` times along the first dimension,
`n` times along the second dimension, and so on. `scale` is a real number that re-scales the lattice constant and atom locations.
"""
function generate_sites(lattice::AbstractLattice{D}, repeats::Vararg{Int,D}; scale = 1.0) where {D}
    scale > 0 || error("scale must be a positive real value.")
    return AtomList(
        _generate_sites((lattice_vectors(lattice)...,), (lattice_sites(lattice)...,), repeats...; scale = scale),
    )
end

############ manipulate sites ###############
"""
    offset_axes(sites::AtomList{D, T}, offsets::Vararg{T,D}) where {D, T}
    offset_axes(offsets...)

Offset the `sites` by distance specified by `offsets`.

```jldoctest; setup=:(using BloqadeLattices)
julia> sites = AtomList([(1.0, 2.0), (10.0, 3.0), (1.0, 12.0), (3.0, 5.0)])
4-element AtomList{2, Float64}:
 (1.0, 2.0)
 (10.0, 3.0)
 (1.0, 12.0)
 (3.0, 5.0)

julia> offset_axes(sites, 1.0, 3.0)
4-element AtomList{2, Float64}:
 (2.0, 5.0)
 (11.0, 6.0)
 (2.0, 15.0)
 (4.0, 8.0)
```
"""
function offset_axes(sites::AtomList{D,T}, offsets::Vararg{T,D}) where {D,T}
    @assert length(offsets) == D
    return AtomList(map(x -> ntuple(i -> x[i] + offsets[i], D), sites.atoms))
end

"""
    rescale_axes(sites::AtomList{D, T}, scale::Real) where {D, T}
    rescale_axes(scale)

Rescale the `sites` by a constant `scale`.

```jldoctest; setup=:(using BloqadeLattices)
julia> sites = AtomList([(1.0, 2.0), (10.0, 3.0), (1.0, 12.0), (3.0, 5.0)])
4-element AtomList{2, Float64}:
 (1.0, 2.0)
 (10.0, 3.0)
 (1.0, 12.0)
 (3.0, 5.0)

julia> rescale_axes(sites, 2.0)
4-element AtomList{2, Float64}:
 (2.0, 4.0)
 (20.0, 6.0)
 (2.0, 24.0)
 (6.0, 10.0)
```
"""
function rescale_axes(sites::AtomList{D,T}, scale::Real) where {D,T}
    return AtomList(map(x -> ntuple(i -> x[i] * scale, D), sites.atoms))
end

"""
    random_dropout(sites::AtomList{D, T}, ratio::Real) where {D, T}
    random_dropout(ratio)

Randomly drop out `ratio * number of sites` atoms from `sites`, where `ratio` ∈ [0, 1].
"""
function random_dropout(sites::AtomList{D,T}, ratio::Real) where {D,T}
    (ratio >= 0 && ratio <= 1) || throw(ArgumentError("dropout ratio be in range [0, 1], got `$ratio`."))
    atoms = sample(1:length(sites), round(Int, length(sites) * (1 - ratio)); replace = false)
    return sites[sort!(atoms)]
end

"""
    clip_axes(sites::AtomList{D, T}, bounds::Vararg{Tuple{T,T},D}) where {D, T}
    clip_axes(bounds...)

Remove sites out of `bounds`, where `bounds` is specified by D D-tuples.

```jldoctest; setup=:(using BloqadeLattices)
julia> sites = AtomList([(1.0, 2.0), (10.0, 3.0), (1.0, 12.0), (3.0, 5.0)])
4-element AtomList{2, Float64}:
 (1.0, 2.0)
 (10.0, 3.0)
 (1.0, 12.0)
 (3.0, 5.0)

julia> clip_axes(sites, (-5.0, 5.0), (-5.0, 5.0))
2-element AtomList{2, Float64}:
 (1.0, 2.0)
 (3.0, 5.0)
```
"""
function clip_axes(sites::AtomList{D,T}, bounds::Vararg{Tuple{T,T},D}) where {D,T}
    @assert length(bounds) == D
    @assert all(x -> length(x) == 2, bounds)
    return AtomList(filter(x -> all(i -> bounds[i][1] <= x[i] <= bounds[i][2], 1:D), sites.atoms))
end
clip_axes(args::Vararg{T,D}) where {T,D} = ls -> clip_axes(ls, args...)
offset_axes(args::Vararg{T,D}) where {T,D} = ls -> offset_axes(ls, args...)
random_dropout(probability::Real) = ls -> random_dropout(ls, probability)
rescale_axes(scale::Real) = ls -> rescale_axes(ls, scale)

############ manipulate grid ###############
"""
    MaskedGrid{T}
    MaskedGrid(xs, ys, mask)

Masked square lattice contains 3 fields, the x-coordinates, y-coordinates and a mask.
e.g. `MaskedGrid([0.0, 1.0, 3.0], [0.0, 2.0,6.0], Bool[1 0 0; 0 1 1; 0 1 0])` specifies the following lattice:

         y₁   y₂        y₃
         ↓    ↓         ↓
    x₁ → ●    ⋅         ●
    x₂ → ⋅    ●         ●

    x₃ → ⋅    ●         ⋅
"""
struct MaskedGrid{T}
    xs::Vector{T}
    ys::Vector{T}
    mask::Matrix{Bool}
end

padydim(al::AtomList{1,T}) where {T} = AtomList([(x[1], zero(T)) for x in al.atoms])
padydim(al::AtomList{2,T}) where {T} = al

"""
    make_grid(sites::AtomList; atol=...)

Create a [`MaskedGrid`](@ref) from the sites. It is required by lattice preparation of Rydberg array.
Because the grid will sort the sites by rows, we need `atol` (default value is 10 time sit data precision)
determines up to what level of round off error, two atoms belong to the same row.
"""
function make_grid(sites::AtomList{D,T}; atol = 10 * eps(T)) where {D,T}
    sites = padydim(sites)
    xs = sort!(approximate_unique(getindex.(sites, 1), atol))
    ys = sort!(approximate_unique(getindex.(sites, 2), atol))
    ixs = map(s -> findfirst(==(s[1]), xs), sites)
    iys = map(s -> findfirst(==(s[2]), ys), sites)
    m, n = length(xs), length(ys)
    mask = zeros(Bool, m, n)
    for (ix, iy) in zip(ixs, iys)
        mask[ix, iy] = true
    end
    return MaskedGrid(xs, ys, mask)
end

# return `(uxs, ixs)``, where `uxs` is the unique x-coordinates, `ixs` the mapping from the index in `xs` to the index in `uxs`.
function approximate_unique(xs::AbstractVector{T}, atol) where {T}
    uxs = T[]
    for x in xs
        found = false
        for ux in uxs
            if isapprox(x, ux; atol = atol)
                found = true
                break
            end
        end
        if !found
            push!(uxs, x)
        end
    end
    return uxs
end

"""
    collect_atoms(maskedgrid::MaskedGrid)

Returns an list of atoms in the `maskedgrid` in order.
"""
function collect_atoms(mg::MaskedGrid)
    return AtomList(map(ci -> (mg.xs[ci.I[1]], mg.ys[ci.I[2]]), findall(mg.mask)))
end

# TODO
# pseudo-lattices,
# image/svg output (maybe),


