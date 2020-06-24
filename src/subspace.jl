export Subspace

"""
    Subspace

A `Dict`-like object stores the mapping between subspace and full space.
"""
struct Subspace
    map::OrderedDict{Int, Int} # fullspace_index => subspace_index
    subspace_v::Vector{Int}
end

function Subspace(subspace_v::Vector{Int})
    subspace_v = sort(subspace_v)
    map = OrderedDict{Int, Int}()
    for (subspace_index, fullspace_index) in enumerate(subspace_v)
        map[fullspace_index] = subspace_index
    end
    return Subspace(map, subspace_v)
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
Base.keys(s::Subspace) = keys(s.map)
Base.values(s::Subspace) = values(s.map)
Base.length(s::Subspace) = length(s.subspace_v)
Base.iterate(s::Subspace) = iterate(s.map)
Base.iterate(s::Subspace, st) = iterate(s.map, st)
Base.haskey(s::Subspace, key) = haskey(s.map, key)
