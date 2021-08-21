abstract type AbstractSpace end
struct FullSpace <: AbstractSpace end

"""
    Subspace

A `Dict`-like object stores the mapping between subspace and full space.
"""
struct Subspace{S <: AbstractVector{Int}} <: AbstractSpace
    map::OrderedDict{Int, Int} # fullspace_index => subspace_index
    subspace_v::S
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
    blockade_subspace(graph)

Create a subspace from given graph's maximal independent set.
"""
function blockade_subspace(graph::SimpleGraph)
    if isempty(edges(graph))
        @warn "graph has empty edges, creating a subspace contains the entire fullspace, consider using a full space register."
    end

    cg = complement(graph)
    mis = maximal_cliques(cg)
    n = nv(graph)
    return create_subspace_from_mis(n, mis)
end

"""
    create_subspace_from_mis(n::Int, mis::AbstractVector)

Create `Subspace` from given list of maximal cliques/maximal independent set.

# Arguments

- `n`: number of vertices of the graph.
- `mis`: the list of maximal independent set.
"""
function create_subspace_from_mis(n::Int, mis::AbstractVector)
    iterators = ThreadsX.map(mis) do each
        fixed_points = setdiff(1:n, each)
        locs = sort!(fixed_points)
        # copied from bsubspace to use runtime length
        masks, shift_len = group_shift(locs)
        len = 1 << (n - length(locs))
        BitSubspace(n, len, length(masks), masks, shift_len)
    end

    # NOTE: ThreadsX doesn't support auto init when using union as op
    # need the following PR merged
    # https://github.com/JuliaFolds/InitialValues.jl/pull/60
    subspace_v = ThreadsX.reduce(union, iterators; init=Int[])
    return Subspace(subspace_v)
end

"""
    blockade_subspace(atoms[, radius=1.0])

Create a blockade approximation subspace from given atom positions and radius.
"""
function blockade_subspace(atoms::AbstractVector{<:RydAtom}, radius::AbstractFloat=1.0)
    return blockade_subspace(unit_disk_graph(atoms, radius))
end

Base.getindex(s::Subspace, key::Int) = s.map[key]
Base.getindex(s::Subspace, key::BitStr) = s.map[buffer(key)]
Base.keys(s::Subspace) = keys(s.map)
Base.values(s::Subspace) = values(s.map)
Base.length(s::Subspace) = length(s.subspace_v)
Base.iterate(s::Subspace) = iterate(s.map)
Base.iterate(s::Subspace, st) = iterate(s.map, st)
Base.haskey(s::Subspace, key) = haskey(s.map, key)
Base.copy(s::Subspace) = Subspace(copy(s.map), copy(s.subspace_v))
Base.vec(s::Subspace) = s.subspace_v
Base.:(==)(x::Subspace, y::Subspace) = (x.map == y.map) && (x.subspace_v == y.subspace_v)
