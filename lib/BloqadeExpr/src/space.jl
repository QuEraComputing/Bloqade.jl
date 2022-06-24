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
    Subspace{S <: AbstractVector} <: AbstractSpace

A `Dict`-like object stores the mapping between subspace and full space.
"""
struct Subspace{T, S<:AbstractVector{T}} <: AbstractSpace
    nqubits::Int
    map::Dict{T,T} # fullspace_index => subspace_index
    subspace_v::S

    function Subspace(nqubits::Int, map::Dict{T,T}, subspace_v) where T<:Integer
        maximum(subspace_v) ≤ one(T) << nqubits || throw(ArgumentError("subspace index is too large"))
        return new{T,typeof(subspace_v)}(nqubits, map, subspace_v)
    end
end

"""
    Subspace(nqubits::Int, subspace_v::AbstractVector)

Create a Subspace from given list of subspace indices in the corresponding full space.
"""
function Subspace(nqubits::Int, subspace_v::AbstractVector{T}) where T<:Integer
    subspace_v = sort(subspace_v)
    map = Dict{T,T}()
    for (subspace_index, fullspace_index) in enumerate(subspace_v)
        map[fullspace_index] = subspace_index
    end
    return Subspace(nqubits, map, subspace_v)
end

function Base.show(io::IO, ::MIME"text/plain", s::Subspace)
    print(io, s.nqubits, "-qubits ", length(s.subspace_v), "-elements ")
    summary(io, s)
    println(io, ":")

    mid_indent = 2
    tab(n) = " "^n

    N = length(s.subspace_v)
    lcol_width = ndigits(N) + mid_indent
    rcol_title = "fullspace"

    function print_line(sub_idx, ful_idx)
        print(io, tab(ndigits(N) - ndigits(sub_idx) + 1))
        printstyled(io, sub_idx; color = :light_black)
        # NOTE: right col also need a tab
        print(io, '│', tab(1))
        return printstyled(io, ful_idx)
    end

    println(io, '─'^(ndigits(N) + 1), "┬", '─'^(ndigits(maximum(s.subspace_v)) + 1))

    if get(io, :limit, false) && N > 10 # compact print
        print_line(1, s.subspace_v[1])
        println(io)
        print_line(2, s.subspace_v[2])
        println(io)
        print_line(3, s.subspace_v[3])
        println(io)
        print(io, tab(ndigits(N) ÷ 2 + 1))
        printstyled(io, tab(1), '⋮'; color = :light_black)
        println(io, '│', tab(1), '⋮')
        print_line(N - 2, s.subspace_v[end-2])
        println(io)
        print_line(N - 1, s.subspace_v[end-1])
        println(io)
        print_line(N, s.subspace_v[end])
    else
        for (sub_idx, ful_idx) in enumerate(s.subspace_v)
            print_line(sub_idx, ful_idx)
            if sub_idx != N
                println(io)
            end
        end
    end
end

Base.getindex(s::Subspace, key::Integer) = s.map[key]
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
    return (x.nqubits == y.nqubits) && (x.map == y.map) && (x.subspace_v == y.subspace_v)
end
