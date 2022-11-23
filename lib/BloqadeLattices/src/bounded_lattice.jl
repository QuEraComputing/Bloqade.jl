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
    parallelepiped_region(lattice::AbstractLattice{D},M::Vararg{NTuple{D,Int},D};pbc::Bool=false)

Create a `BoundedLattice` given an existing lattice and tuples defining a parallelogram/paralelepiped
/line segment defined by vectors that are integer multiples of the lattice vectors in `lattice`.

Periodic Boundary Conditions can be enable/disabled via `pbc`.

Tuples must be the same length and quantity as the dimensions of the lattice argument.

```jldoctest; setup=:(using BloqadeLattices)
julia>  parallelepiped_region(SquareLattice(),(2,0),(0,2);pbc=true)
BoundedLattice{SquareLattice, Parallelepiped{2, Float64}}(SquareLattice(), Parallelepiped{2, Float64}([2.0 0.0; 0.0 2.0], [0.5 0.0; 0.0 0.5]), [(0.0, 0.0), (0.0, 1.0), (1.0, 0.0), (1.0, 1.0)], true)

julia> parallelepiped_region(KagomeLattice(),(2,2),(-2,2))
BoundedLattice{KagomeLattice, Parallelepiped{2, Float64}}(KagomeLattice(), Parallelepiped{2, Float64}([3.0 -1.0; 1.7320508075688772 1.7320508075688772], [0.25 0.14433756729740643; -0.25 0.43301270189221935]), [(-0.75, 1.299038105676658), (-0.5, 0.8660254037844386), (-0.25, 0.4330127018922193), (-0.25, 1.299038105676658), (0.0, 0.0), (0.0, 1.7320508075688772), (0.25, 0.4330127018922193), (0.25, 1.299038105676658), (0.25, 2.1650635094610964), (0.5, 0.8660254037844386)  …  (1.25, 1.299038105676658), (1.25, 2.1650635094610964), (1.5, 0.8660254037844386), (1.5, 2.598076211353316), (1.75, 1.299038105676658), (1.75, 2.1650635094610964), (1.75, 3.031088913245535), (2.0, 1.7320508075688772), (2.25, 1.299038105676658), (2.25, 2.1650635094610964)], false)
```
"""
function parallelepiped_region(lattice::AbstractLattice{D},M::Vararg{NTuple{D,Int},D};pbc::Bool=false) where D
    lat_vecs = lattice_vectors(lattice)
    T = eltype(lat_vecs[1])
    bounds =  zeros(T,D,D)
    
    for i in 1:D
        for j in 1:D
            bounds[:,i] .+= M[i][j] .* lat_vecs[j]
        end
    end
    region = Parallelepiped(bounds)
    
    return BoundedLattice(lattice,region,pbc)

end

"""
    dimension(lattice::BoundedLattice{L,C})

Returns the dimensions of the `BoundedLattice` (e.g.: `2` for 2D, `3` for 3D)

```jldoctest; setup=:(using BloqadeLattices)
julia> bl = parallelepiped_region(ChainLattice(),(4,);pbc=true) # create a 1D BoundedLattice
BoundedLattice{ChainLattice, Parallelepiped{1, Float64}}(ChainLattice(), Parallelepiped{1, Float64}([4.0;;], [0.25;;]), [(0.0,), (1.0,), (2.0,), (3.0,)], true)

julia> dimension(bl)
1

julia> bl = parallelepiped_region(SquareLattice(),(3,0),(0,2);) # create a 2D BoundedLattice
BoundedLattice{SquareLattice, Parallelepiped{2, Float64}}(SquareLattice(), Parallelepiped{2, Float64}([3.0 0.0; 0.0 2.0], [0.3333333333333333 0.0; 0.0 0.5]), [(0.0, 0.0), (0.0, 1.0), (1.0, 0.0), (1.0, 1.0), (2.0, 0.0), (2.0, 1.0)], false)

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
BoundedLattice{SquareLattice, Parallelepiped{2, Float64}}(SquareLattice(), Parallelepiped{2, Float64}([3.0 0.0; 0.0 2.0], [0.3333333333333333 0.0; 0.0 0.5]), [(0.0, 0.0), (0.0, 1.0), (1.0, 0.0), (1.0, 1.0), (2.0, 0.0), (2.0, 1.0)], false)

julia> lattice_vectors(bl) # lattice vectors used in Bravais Lattice definition of underlying SquareLattice
((1.0, 0.0), (0.0, 1.0))
```
"""
lattice_vectors(lattice::BoundedLattice{L,C}) where {L,C} = lattice_vectors(lattice.lattice)

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
BoundedLattice{SquareLattice, Parallelepiped{2, Float64}}(SquareLattice(), Parallelepiped{2, Float64}([1.0 0.0; 0.0 1.0], [1.0 0.0; 0.0 1.0]), [(0.0, 0.0)], false)

julia> distance(bl, (0.1, 0.1), (0.5, 1.1)) # distance between two points
1.077032961426901

julia> bl_pbc = parallelepiped_region(SquareLattice(), (1,0),(0,1);pbc=true) # Define 2D BoundedLattice with Periodic Boundary Condition
BoundedLattice{SquareLattice, Parallelepiped{2, Float64}}(SquareLattice(), Parallelepiped{2, Float64}([1.0 0.0; 0.0 1.0], [1.0 0.0; 0.0 1.0]), [(0.0, 0.0)], true)

julia> distance(bl_pbc, (0.1, 0.1), (0.5, 1.1)) # distance with periodic boundary condition enabled
0.4
```
"""
distance(lattice::BoundedLattice,x,y) = lattice.pbc ? distance(lattice.region,x,y) : distance(x,y)
