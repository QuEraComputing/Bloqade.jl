# abstraction for a single tile of an infinite lattice
"""
    struct BoundedLattice{L<:AbstractLattice, R<:AbstractRegion}
    
Defines a lattice bounded by a region with the option for periodic boundary conditions.

# Fields
- `lattice <: AbstractLattice`: Lattice to be bounded.
- `region <: AbstractRegion`: Region that bounds the underlying `lattice`.
- `site_positions::AtomList`: Positions of the atoms inside `region`.
- `pbc::Bool`: Enable/Disable behavior for Periodic Boundary Conditions.

"""
struct BoundedLattice{L<:AbstractLattice,R<:AbstractRegion} 
    lattice::L
    region::R
    site_positions::AtomList
    pbc::Bool
    function BoundedLattice(lattice::L,region::R,site_positions,pbc::Bool) where {L<:AbstractLattice,R<:AbstractRegion} 
        sort!(site_positions)

        return new{L,R}(lattice,region,AtomList(site_positions),pbc)
    end
end

"""
    BoundedLattice(lattice::AbstractLattice{D},region::AbstractRegion{D},pbc::Bool=false)

Creates a [`BoundedLattice`](@ref) instance when provided with the underlying `lattice` and `region` to 
bound on the lattice, with the option to enable Periodic Boundary Conditions.

See also: [`parallelepiped_region`](@ref)
"""
function BoundedLattice(lattice::AbstractLattice{D},region::AbstractRegion{D},pbc::Bool=false) where D
    site_positions = generate_sites_in_region(lattice,region)
    return BoundedLattice(lattice,region,site_positions,pbc)
end

"""
    parallelepiped_region(lattice::AbstractLattice{D},M::Vararg{NTuple{D,Int},D};pbc::Bool=false;scale::Real=1)

Create a `BoundedLattice` given an existing lattice and tuples defining a parallelogram/paralelepiped
/line segment defined by vectors that are integer multiples of the lattice vectors in `lattice`.

Periodic Boundary Conditions can be enable/disabled via `pbc`.
Tuples must be the same length and quantity as the dimensions of the lattice argument.

```jldoctest; setup=:(using BloqadeLattices)
julia> parallelepiped_region(SquareLattice(),(2,0),(0,2);pbc=true);

julia> parallelepiped_region(KagomeLattice(),(2,2),(-2,2));
```
"""
function parallelepiped_region(lattice::AbstractLattice{D},M::Vararg{NTuple{D,Int},D};pbc::Bool=false,scale::Real=1.0) where D
    scale > 0 || error("scale must be a positive real value.")

    lat_vecs = lattice_vectors(lattice)
    lat_sites = lattice_sites(lattice)
    T = promote_type(eltype(lat_vecs[1]),Int,typeof(scale))
    scaled_bounds =  zeros(T,D,D)
    scaled_lat_vecs = Tuple(scale .* lat_vec for lat_vec in lat_vecs)
    scaled_lat_sites = Tuple(scale .* lat_site for lat_site in lat_sites)

    for i in 1:D
        for j in 1:D
            scaled_bounds[:,i] .+= M[i][j] .* scaled_lat_vecs[j]
        end
    end

    scaled_region = Parallelepiped(scaled_bounds)
    scaled_lattice = GeneralLattice(scaled_lat_vecs,scaled_lat_sites)
    
    return BoundedLattice(scaled_lattice,scaled_region,pbc)

end

"""
    dimension(lattice::BoundedLattice{L,C})

Returns the dimensions of the `BoundedLattice` (e.g.: `2` for 2D, `3` for 3D)

```jldoctest; setup=:(using BloqadeLattices)
julia> bl = parallelepiped_region(ChainLattice(),(4,);pbc=true) # create a 1D BoundedLattice

julia> dimension(bl)
1

julia> bl = parallelepiped_region(SquareLattice(),(3,0),(0,2);) # create a 2D BoundedLattice

julia> dimension(bl)
2
```
"""
dimension(lattice::BoundedLattice{L,C}) where {L,C} = dimension(lattice.lattice)

"""
    lattice_vectors(lattice::BoundedLattice{L,C})

Returns the underlying Bravais lattice vectors of the `BoundedLattice`

```jldoctest; setup=:(using BloqadeLattices)
julia> bl = parallelepiped_region(SquareLattice(),(3,0),(0,2);) # create a 2D BoundedLattice

julia> lattice_vectors(bl) # lattice vectors used in Bravais Lattice definition of underlying SquareLattice
((1.0, 0.0), (0.0, 1.0))
```
"""
lattice_vectors(lattice::BoundedLattice{L,C}) where {L,C} = lattice_vectors(lattice.lattice)

"""
    lattice_sites(lattice::BoundedLattice{L,C})

Returns the underlying site vectors that define the unit-cell of the `BoundedLattice`

```jldoctest; setup=:(using BloqadeLattices)
julia> bl = parallelepiped_region(SquareLattice(),(3,0),(0,2);) # create a 2D BoundedLattice

julia> lattice_sites(bl) # lattice vectors used in Bravais Lattice definition of underlying SquareLattice
((0.0, 0.0), )
```
"""
lattice_sites(lattice::BoundedLattice{L,C}) where {L,C} = lattice_sites(lattice.lattice)

Base.length(lattice::BoundedLattice{L,C}) where {L,C} = length(lattice.site_positions)

# what behavior to produce when position is not found?
function get_position_index(pos,lattice::BoundedLattice{L,C}) where {L,C}
    j = searchsortedlast(Tuple(pos),lattice.site_positions)
    !(pos ≈ lattice.site_positions[j]) && error("$pos not contained in bounded lattice.")
    return j
end

"""
    distance(lattice::BoundedLattice,x,y)

Returns the distance between two points in the [`BoundedLattice`](@ref). 

Points `x` and `y` can be any iterable and must have the same dimensions as the [`BoundedLattice`](@ref) 
(ex: `(x,y)` for a 2D lattice, `(x,y,z)` for a 3D lattice).

If the Periodic Boundary Condition option has been set to `true` for the [`BoundedLattice`](@ref),
the smallest distance between points (modulo the region) is returned, otherwise the standard Euclidean
metric is used.

```jldoctest; setup=:(using BloqadeLattices)
julia> bl = parallelepiped_region(SquareLattice(), (1,0),(0,1);) # Define 2D BoundedLattice

julia> distance(bl, (0.1, 0.1), (0.5, 1.1)) # distance between two points
1.077032961426901

julia> bl_pbc = parallelepiped_region(SquareLattice(), (1,0),(0,1);pbc=true) # Define 2D BoundedLattice with Periodic Boundary Condition

julia> distance(bl_pbc, (0.1, 0.1), (0.5, 1.1)) # distance with periodic boundary condition enabled
0.4
```
"""
distance(lattice::BoundedLattice,x,y) = lattice.pbc ? distance(lattice.region,x,y) : distance(x,y)
