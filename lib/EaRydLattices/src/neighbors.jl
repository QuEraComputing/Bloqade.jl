"""
    make_kdtree(locations::AbstractVector{NTuple{D,T}}) where {T, D}

Returns a [`KDTree`](@ref) instance from input `locations`.
"""
function make_kdtree(locations::AbstractVector{NTuple{D,T}}) where {T, D}
    data = zeros(T, D, length(locations))
    for i=1:length(locations)
        data[:,i] .= locations[i]
    end
    return KDTree(data)
end

"""
    DistanceGroup

The vertices grouped by distances. One can use `distancegroup[n]` to get `n`-th nearest neighbors.
"""
struct DistanceGroup
    siteindices::Vector{Int}
    ptr::Vector{Int}
end
function Base.getindex(dg::DistanceGroup, k::Int)
    if length(dg.ptr) < k+2
        return Int[]
    else
        return dg.siteindices[dg.ptr[k+1]:dg.ptr[k+2]-1]
    end
end

# Obtain `K` vertices closest to the site specified by `siteindex` from the KDTree `tree`.
# It returns a 2-tuple of vertex indices and distances.
function nearest(tree::KDTree, siteindex::Int, K::Int)
    siteindices, distances = knn(tree, tree.data[findfirst(==(siteindex), tree.indices)], K)
    perm = sortperm(distances)
    return siteindices[perm], distances[perm]
end

# Group the `siteindices` by their `distances`. Returns a `DistanceGroup` instance that can be used to index n-th nearest neighbors.
function group_by_distances(siteindices, distances::AbstractVector{T}, atol) where T
    dpre = distances[1]
    ptr = Int[1]
    for i=2:length(siteindices)
        di = distances[i]
        if !isapprox(di, dpre; atol=atol)
            push!(ptr, i)  # endpoint
            dpre = di
        end
    end
    push!(ptr, length(siteindices)+1)  # endpoint
    return DistanceGroup(siteindices, ptr)
end

"""
    grouped_nearest(tree::KDTree, siteindex::Int, K::Int; atol=1e-8)

Find the `K` closest vertices to `siteindex`, and group them by distance (distances different by a number smaller than 1e-8 are treated as the same).
Returns a [`DistanceGroup`](@ref) instance.

```jldoctest; setup=:(using EaRydLattices)
julia> locations = generate_sites(HoneycombLattice(), 5, 5);

julia> tree = make_kdtree(locations)
NearestNeighbors.KDTree{StaticArrays.SVector{2, Float64}, Distances.Euclidean, Float64}
  Number of points: 50
  Dimensions: 2
  Metric: Distances.Euclidean(0.0)
  Reordered: true

julia> gn = grouped_nearest(tree, 23, 20)
DistanceGroup([23, 14, 22, 24, 15, 13, 21, 25, 33, 31, 12, 16, 32, 4, 6, 34, 26, 17, 5, 41], [1, 2, 5, 11, 14, 18, 21])

julia> gn[0]  # the 0-th nearest neighbor is defined by vertex itself
1-element Vector{Int64}:
 23

julia> gn[1]  # nearest neighbors
3-element Vector{Int64}:
 14
 22
 24

julia> gn[2]  # second nearest neighbors
6-element Vector{Int64}:
 15
 13
 21
 25
 33
 31
```
"""
function grouped_nearest(tree::KDTree, siteindex::Int, K::Int; atol=1e-8)
    group_by_distances(nearest(tree, siteindex, K)..., atol)
end