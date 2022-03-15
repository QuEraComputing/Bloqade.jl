"""
    AbstractSpace

Abstract type for spaces.
"""
abstract type AbstractSpace end

"""
    FullSpace <: AbstractSpace

A trait for the full space.
"""
struct FullSpace <: AbstractSpace end

"""
    fullspace

A constant for the [`FullSpace`](@ref).
"""
const fullspace = FullSpace()

"""
    Subspace{S <: AbstractVector{Int}} <: AbstractSpace

A `Dict`-like object stores the mapping between subspace and full space.
"""
struct Subspace{S <: AbstractVector{Int}} <: AbstractSpace
    nqubits::Int
    map::Dict{Int, Int} # fullspace_index => subspace_index
    subspace_v::S
end

"""
    Subspace(nqubits::Int, subspace_v::AbstractVector{Int})

Create a Subspace from given list of subspace indices in the corresponding full space.
"""
function Subspace(nqubits::Int, subspace_v::AbstractVector{Int})
    subspace_v = sort(subspace_v)
    map = Dict{Int, Int}()
    for (subspace_index, fullspace_index) in enumerate(subspace_v)
        map[fullspace_index] = subspace_index
    end
    return Subspace(nqubits, map, subspace_v)
end

function Base.show(io::IO, ::MIME"text/plain", s::Subspace{S}) where S
    print(io, s.nqubits, "-qubits ", length(s.subspace_v), "-elements ")
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
        printstyled(io, tab(1), '⋮'; color=:light_black)
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

Base.getindex(s::Subspace, key::Int) = s.map[key]
Base.getindex(s::Subspace, key::BitStr) = s.map[buffer(key)]
Base.keys(s::Subspace) = keys(s.map)
Base.values(s::Subspace) = values(s.map)
Base.length(s::Subspace) = length(s.subspace_v)
Base.iterate(s::Subspace) = iterate(s.map)
Base.iterate(s::Subspace, st) = iterate(s.map, st)
Base.haskey(s::Subspace, key) = haskey(s.map, key)
Base.copy(s::Subspace) = Subspace(s.nqubits, copy(s.map), copy(s.subspace_v))
Base.vec(s::Subspace) = s.subspace_v
Base.to_index(ss::Subspace) = ss.subspace_v .+ 1

function Base.:(==)(x::Subspace, y::Subspace)
    return (x.nqubits == y.nqubits) && (x.map == y.map) &&
        (x.subspace_v == y.subspace_v)
end

# TODO: move these to MIS
# """
#     independent_set_subspace(graph)

# Create a subspace from given graph's maximal independent set.
# """
# function independent_set_subspace(graph::SimpleGraph)
#     if isempty(edges(graph))
#         @warn "graph has empty edges, creating a subspace contains the entire fullspace, consider using a full space register."
#     end

#     cg = complement(graph)
#     mis = maximal_cliques(cg)
#     n = nv(graph)
#     return create_subspace_from_mis(n, mis)
# end

# """
#     create_subspace_from_mis(n::Int, mis::AbstractVector)

# Create `Subspace` from given list of maximal cliques/maximal independent set.

# # Arguments

# - `n`: number of vertices of the graph.
# - `mis`: the list of maximal independent set.
# """
# function create_subspace_from_mis(n::Int, mis::AbstractVector)
#     iterators = ThreadsX.map(mis) do each
#         fixed_points = setdiff(1:n, each)
#         locs = sort!(fixed_points)
#         # copied from bsubspace to use runtime length
#         masks, shift_len = group_shift(locs)
#         len = 1 << (n - length(locs))
#         BitSubspace(n, len, length(masks), masks, shift_len)
#     end

#     # NOTE: ThreadsX doesn't support auto init when using union as op
#     # need the following PR merged
#     # https://github.com/JuliaFolds/InitialValues.jl/pull/60
#     subspace_v = ThreadsX.reduce(union, iterators; init=OnInit(Vector{Int}))
#     return Subspace(n, subspace_v)
# end

# """
#     blockade_subspace(atoms[, radius=1.0])

# Create a blockade approximation subspace from given atom positions and radius.
# """
# function blockade_subspace(atoms::AbstractVector{<:RydAtom}, radius::AbstractFloat=1.0)
#     return independent_set_subspace(unit_disk_graph(atoms, radius))
# end

# function blockade_subspace(atoms::AbstractVector{<:Tuple}, radius::AbstractFloat=1.0)
#     blockade_subspace(RydAtom.(atoms), radius)
# end
