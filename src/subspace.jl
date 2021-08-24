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

function Base.show(io::IO, ::MIME"text/plain", s::Subspace{S}) where S
    print(io, length(s.subspace_v), "-elements ")
    summary(io, s)
    println(io, ":")

    mid_indent = 2
    tab(n) = " "^n

    N = length(s.subspace_v)
    lcol_width = ndigits(N) + mid_indent
    rcol_title = "fullspace"

    function print_line(sub_idx, ful_idx)
        print(io, tab(ndigits(N)-ndigits(sub_idx)+1))
        printstyled(io, sub_idx; color=:light_black)
        # NOTE: right col also need a tab
        print(io, '│', tab(1))
        printstyled(io, ful_idx)
    end

    println(io, '─'^(ndigits(N)+1), "┬", '─'^(ndigits(maximum(s.subspace_v))+1))

    if get(io, :limit, false) && N > 10 # compact print
        print_line(1, s.subspace_v[1])
        println(io)
        print_line(2, s.subspace_v[2])
        println(io)
        print_line(3, s.subspace_v[3])
        println(io)
        print(io, tab(ndigits(N)÷2+1))
        printstyled(io, '⋮'; color=:light_black)
        println(io, '│', tab(1), '⋮')
        print_line(N-2, s.subspace_v[end-2])
        println(io)
        print_line(N-1, s.subspace_v[end-1])
        println(io)
        print_line(N  , s.subspace_v[end  ])
    else
        for (sub_idx, ful_idx) in enumerate(s.subspace_v)
            print_line(sub_idx, ful_idx)
            if sub_idx != N
                println(io)
            end
        end
    end
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
    subspace_v = ThreadsX.reduce(union, iterators; init=OnInit(Vector{Int}))
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
