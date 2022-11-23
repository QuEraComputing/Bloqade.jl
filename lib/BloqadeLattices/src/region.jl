"""
    AbstractRegion{D}

Supertype for all `D` dimensional regions used to define bounds on lattices.

# Implementation

The following should be overriden:

- `Base.in`: Returns `Bool` on whether a point exists in the region
- `Base.mod`: Maps points outside the region back into the region
- `distance`: Calculates distance between points with periodic boundary conditions enabled
"""
abstract type AbstractRegion{D} end

Base.broadcastable(x::AbstractRegion) = Ref(x)

function generate_neighboring_sites(site::NTuple{N,Int},n_sites::Int) where N
   
    # site = (n[1],n[2],...n[D],i)
    neighboring_sites = NTuple{N,Int}[]
    for i in 1:N-1 # add sites from neighboring unit cells
        new_site = [site...,]
        
        for k in 1:n_sites 

            new_site[end] = k
            new_site[i] += 1
            push!(neighboring_sites,Tuple(new_site))
            new_site .= site

            new_site[end] = k
            new_site[i] -= 1
            push!(neighboring_sites,Tuple(new_site))
            new_site .= site            
            
        end
    end

    return neighboring_sites
end

"""
    generate_sites_in_region(lattice::AbstractLattice{D}, region::AbstractRegion{D})

Generates sites from the `lattice` that are present in the `region`.

```jldoctest; setup=:(using BloqadeLattices)
julia> generate_sites_in_region(ChainLattice(), Parallelepiped(4.0))
4-element Vector{Tuple{Float64}}:
 (0.0,)
 (1.0,)
 (2.0,)
 (3.0,)

julia> bounds = zeros(2,2)
2×2 Matrix{Float64}:
 0.0  0.0
 0.0  0.0

julia> bounds[:,1] .= (0.0, 2.0); bounds[:,2] .= (2.0, 0.0);

julia> generate_sites_in_region(SquareLattice(), Parallelepiped(bounds))
4-element Vector{Tuple{Float64, Float64}}:
 (0.0, 0.0)
 (1.0, 0.0)
 (0.0, 1.0)
 (1.0, 1.0)
```
"""
function generate_sites_in_region(lattice::AbstractLattice{D}, region::AbstractRegion{D}) where D
    zeros(D) ∉ region && error("bounding region must contain origin")

    lat_vecs = lattice_vectors(lattice)
    lat_sites = lattice_sites(lattice)
    T = eltype(first(lat_vecs))

    origin = Tuple(vcat(zeros(Int,D),[1]))
    stack = NTuple{D+1,Int}[origin]
    visited_sites = Set{NTuple{D+1,Int}}([origin])
    site_positions = NTuple{D,T}[Tuple(zeros(D))]

    while !isempty(stack)
        site = pop!(stack)
        
        for nn_site in generate_neighboring_sites(site,length(lat_sites))
            nn_site_position = sum([lv...,] .* n for (lv,n) in zip(lat_vecs,nn_site)) .+ lat_sites[nn_site[end]]
            if nn_site ∉ visited_sites && nn_site_position ∈ region
                push!(stack, nn_site)
                push!(visited_sites, nn_site)
                push!(site_positions,Tuple(nn_site_position))
            end
        end
    end

    return site_positions

end

"""
    struct Parallelepiped{D, T} <: AbstractRegion{D}
    
Region that is a Parallelogram/Parallelepiped

See also [`Parallelepiped(vecs)`](@ref), [`Parallelepiped(vecs::T) where {T<:Real}`](@ref)

# Fields
- `vecs::Matrix{T}`: Matrix with column vectors defining Parallelogram/Parallelepiped 
- `vecs_inv::Matrix{T}`: Inverse of `vecs`.
"""
struct Parallelepiped{D, T} <: AbstractRegion{D}
    # abstract matrix dimensions already guaranteed to be 2D
    vecs::Matrix{T}
    vecs_inv::Matrix{T}
end

"""
    Parallelepiped(vecs)
    Parallelepiped(vecs::T) where {T<:Real}

Define a region (either a line segment, parallelogram, or parallelepiped depending
on the dimensions of `vecs`) using a single value or column vectors in a matrix that can be
used to create a [`BoundedLattice`](@ref).

```jldoctest; setup=:(using BloqadeLattices)
julia> Parallelepiped(2.0) # can bound a 1D lattice
Parallelepiped{1, Float64}([2.0;;], [0.5;;])

julia> bounds = zeros((2,2)); # Create 2x2 matrix to store vectors defining 2D region

julia> bounds[:,1] .= (3,3); bounds[:,2] .= (4,0); # Column Vectors define the Parallelogram

julia> Parallelepiped(bounds)
Parallelepiped{2, Float64}([3.0 4.0; 3.0 0.0], [0.0 0.3333333333333333; 0.25 -0.25])
```
"""
function Parallelepiped(vecs)
    D = size(vecs,1)
    vecs_inv = inv(vecs)
    # return Parallelpiped
    return Parallelepiped{D,eltype(vecs_inv)}(vecs,vecs_inv)
end

# Handle the 1D case
function Parallelepiped(vecs::T) where {T<:Real}
    vecs = fill(vecs,1,1)
    vecs_inv = inv(vecs)
    # return Parallelpiped
    return Parallelepiped{1,eltype(vecs_inv)}(vecs,vecs_inv)
end



# check if a point is in the region

approx_in_range(x) = (zero(x) ≤ x < one(x) || isapprox(x,zero(x),atol=eps(typeof(x)))) && !isapprox(x,one(x),atol=eps(typeof(x)))
"""
    Base.in(x,region::Parallelepiped{D,T}) where {D,T} = all(in_range.(region.vecs_inv * [x...,]))

Given a point `x`, check if exists in the `region`.

```jldoctest; setup=:(using BloqadeLattices)
julia> bounds = zeros(2,2);

julia> bounds[:,1] .= (3,3); bounds[:,2] .= (4,0);

julia> p = Parallelepiped(bounds)
Parallelepiped{2, Float64}([3.0 4.0; 3.0 0.0], [0.0 0.3333333333333333; 0.25 -0.25])

julia> (0.0, 0.0) ∈ p
true

julia> (5.0, 1.0) ∈ p
false
```
"""
Base.in(x,region::Parallelepiped{D,T}) where {D,T} = all(approx_in_range.(region.vecs_inv * [x...,]))

"""
    Base.mod(x,region::Parallelepiped{D,T})
    Base.mod(x,region::Parallelepiped{1,T})

Given a point `x`, enforce periodic boundary conditions by having points that fall outside of
the `region` map to ones on the inside.

```jldoctest; setup=:(using BloqadeLattices)
julia> p = Parallelepiped(1.0)

julia> mod(1.0, p) # 1.0 falls outside of region, mapped to opposite end
0.0

julia> bounds = zeros(2,2); bounds[:,1] .= (0.0, 2.0); bounds[:,2] .= (2.0, 0.0); # define square region

julia> mod((1.0, 1.0), Parallelepiped(bounds)) # point is in bounds, so nothing happens
(1.0, 1.0)

julia> mod((3.0, 3.0), Parallelepiped(bounds)) # point is out out of bounds, mapped back inside
(1.0, 1.0)
```
"""
# Enforce periodic boundary conditions by having points that fall outside of the 
# parallelogram map to ones on the inside
Base.mod(x,region::Parallelepiped{D,T}) where {D,T} = typeof(x)(region.vecs * mod.((region.vecs_inv * [x...,]),1.0))
Base.mod(x,region::Parallelepiped{1,T}) where {T} =  mod.(x,only(region.vecs))

"""
    distance(x, y)
    distance(region::Parallelepiped{D,T},x,y)

Distance between two points. 

If just points `x` and `y` are provided, the Euclidean distance is calculated. 

If a `region` is provided, then it is automatically assumed periodic boundary conditions are enabled
and the smallest possible distance between the two points is returned.

```jldoctest; setup=:(using BloqadeLattices)
julia> distance((0.0, 0.0), (1.0, 1.0))
1.4142135623730951

julia> bounds = zeros(2,2); bounds[:,1] .= (3, 0); bounds[:,2] .= (0, 3);

julia> distance(Parallelepiped(bounds), (0.5, 0.5), (2.5, 2.5))
1.4142135623730951
```
"""
distance(x, y) = sqrt(mapreduce(x -> x^2, +, x .- y))

function distance(region::Parallelepiped{D,T},x,y) where {D,T}
    x = mod([x...,],region)
    y = mod([y...,],region)
    dist = distance(x,y)

    for a in 1:D
        shift = region.vecs[:,a]
        for b in a:D
            shift .+= (b > a ? region.vecs[:,b] : 0)
            shift_y = y .+ shift
            shift_x = x .+ shift
            dist = min(
                distance(x,shift_y),
                distance(shift_x,y),
                dist
            )
        end
    end

    return dist
end