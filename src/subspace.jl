export Subspace

struct Subspace
    map::Dict{Int, Int} # fullspace_index => subspace_index
end

function Subspace(subspace_v::Vector{Int})
    subspace_v = sort(subspace_v)
    map = Dict{Int, Int}()
    for (subspace_index, fullspace_index) in enumerate(subspace_v)
        map[fullspace_index] = subspace_index
    end
    return Subspace(map)
end

"""
    Subspace(graph)

Create a subspace from given graph's maximal independent set.
"""
function Subspace(graph::SimpleGraph)
    cg = complement(graph)
    mis = maximal_cliques(cg)
    n = nv(graph)
    it = map(mis) do each
        fixed_points = setdiff(1:n, each)
        itercontrol(n, fixed_points, zero(fixed_points))
    end
    return Subspace(unique(Iterators.flatten(it)))
end

Base.getindex(s::Subspace, key::Int) = s.map[key]
