using NearestNeighbors

export make_kdtree, grouped_nearest

function make_kdtree(locations::AbstractVector{NTuple{D,T}}) where {T, D}
    data = zeros(T, D, length(locations))
    for i=1:length(locations)
        data[:,i] .= locations[i]
    end
    return KDTree(data)
end

struct DistanceGroup
    nodes::Vector{Int}
    ptr::Vector{Int}
end
function Base.getindex(dg::DistanceGroup, k::Int)
    if length(dg.ptr) < k+2
        return Int[]
    else
        return dg.nodes[dg.ptr[k+1]:dg.ptr[k+2]-1]
    end
end

function nearest(tree::KDTree, siteindex::Int, K::Int)
    nodes, distances = knn(tree, tree.data[findfirst(==(siteindex), tree.indices)], K)
    perm = sortperm(distances)
    return nodes[perm], distances[perm]
end

function group_by_distances(nodes, distances::AbstractVector{T}, atol) where T
    dpre = distances[1]
    ptr = Int[1]
    for i=2:length(nodes)
        di = distances[i]
        if !isapprox(di, dpre; atol=atol)
            push!(ptr, i)  # endpoint
            dpre = di
        end
    end
    push!(ptr, length(nodes)+1)  # endpoint
    return DistanceGroup(nodes, ptr)
end

function grouped_nearest(tree::KDTree, siteindex::Int, K::Int; atol=1e-8)
    group_by_distances(nearest(tree, siteindex, K)..., atol)
end