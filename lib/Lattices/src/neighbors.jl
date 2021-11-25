using NearestNeighbors

export make_kdtree, kneighbors

function make_kdtree(locations::AbstractVector{NTuple{D,T}}) where {T, D}
    data = zeros(T, D, length(locations))
    for i=1:length(locations)
        data[:,i] .= locations[i]
    end
    return KDTree(data)
end

function kneighbors(tree::KDTree{VT,ST,T}, k::Int, siteindex::Int; maxn=length(tree.data), atol=eps(T)*10) where {T,ST,VT}
    # get nearest K distances
    K = min(maxn, 6*k+1, length(tree.data))
    locations, distances = knn(tree, tree.data[findfirst(==(siteindex), tree.indices)], K)
    perm = sortperm(distances)
    locations, distances = locations[perm], distances[perm]
    dpre = zero(T)
    ptr = Int[1]
    for i=2:K
        di = distances[i]
        if !isapprox(di, dpre; atol=atol)
            push!(ptr, i)  # endpoint
            dpre = di
        end
    end
    push!(ptr, K+1)  # endpoint
    if length(ptr) <= k
        return Int[]
    else
        return locations[ptr[k+1]:ptr[k+2]-1]
    end
end