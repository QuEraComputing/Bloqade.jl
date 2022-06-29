# copied from BQCESubroutine
# should be removed after we release this package
# mainly because we want to use dynamic sized IterControl here

@inline lmove(b::Integer, mask::Integer, k::Int) = (b & ~mask) << k + (b & mask)

"""
    bmask(::Type{T}) where T <: Integer -> zero(T)
    bmask([T::Type], positions::Int...) -> T
    bmask([T::Type], range::UnitRange{Int}) -> T
Return an integer mask of type `T` where `1` is the position masked according to
`positions` or `range`. Directly use `T` will return an empty mask `0`.
"""
function bmask end

bmask(args...) = bmask(Int, args...)
bmask(::Type{T}) where {T<:Integer} = zero(T)
bmask(::Type{T}, positions::Int...) where {T<:Integer} = bmask(T, positions)

function bmask(::Type{T}, itr) where {T<:Integer}
    isempty(itr) && return 0
    ret = zero(T)
    for b in itr
        ret += one(T) << (b - 1)
    end
    return ret
end

@inline function bmask(::Type{T}, range::UnitRange{Int})::T where {T<:Integer}
    return ((one(T) << (range.stop - range.start + 1)) - one(T)) << (range.start - 1)
end

"""
    anyone(index::Integer, mask::Integer) -> Bool

Return `true` if any masked position of index is 1.

# Example
`true` if any masked positions is 1.
```jldoctest
julia> anyone(0b1011, 0b1001)
true
julia> anyone(0b1011, 0b1100)
true
julia> anyone(0b1011, 0b0100)
false
```
"""
anyone(index::T, mask::T) where {T<:Integer} = (index & mask) != zero(T)

"""
    ismatch(index::Integer, mask::Integer, target::Integer) -> Bool

Return `true` if bits at positions masked by `mask` equal to `1` are equal to `target`.

## Example

```julia
julia> n = 0b11001; mask = 0b10100; target = 0b10000;

julia> ismatch(n, mask, target)
true
```
"""
ismatch(index::T, mask::T, target::T) where {T<:Integer} = (index & mask) == target

@inline function _group_shift(masks::Vector{T}, shift_len::Vector{Int}, k::Int, k_prv::Int) where T<:Integer
    # if current position in the contiguous region
    # since these bits will be moved together with
    # the first one, we don't need to generate a
    # new mask
    if k == k_prv + 1
        shift_len[end] += 1
    else
        # we generate a bit mask where the 1st to k-th bits are 1
        push!(masks, bmask(T, 0:k-1))
        push!(shift_len, 1)
    end
end

@inline function group_shift(::Type{T}, locations) where T<:Integer
    masks = T[]
    shift_len = Int[]
    k_prv = -1
    for k in locations
        _group_shift(masks, shift_len, k, k_prv)
        k_prv = k
    end
    return masks, shift_len
end

struct BitSubspace{L}
    n::Int # number of bits in fullspace
    sz_subspace::L # size of the subspace
    n_shifts::Int # number of shifts
    masks::Vector{L} # shift masks
    shift_len::Vector{Int} # length of each shift
end

function Base.getindex(s::BitSubspace{T}, i::Integer) where T
    index = T(i - 1)
    @inbounds for k in 1:s.n_shifts
        index = lmove(index, s.masks[k], s.shift_len[k])
    end
    return index
end

Base.firstindex(s::BitSubspace) = 1
Base.lastindex(s::BitSubspace) = s.sz_subspace
Base.length(s::BitSubspace) = s.sz_subspace
Base.eltype(::BitSubspace{L}) where L = L

function Base.iterate(s::BitSubspace{T}, st::T = one(T)) where T
    st <= length(s) || return
    return s[st], st + one(T)
end

function Base.show(io::IO, ::MIME"text/plain", s::BitSubspace)
    indent = get(io, :indent, 0)
    println(io, " "^indent, s.sz_subspace, "-element BitSubspace:")
    if s.sz_subspace < 5
        for k in 1:s.sz_subspace
            print(io, " "^(indent + 1), string(s[k]; base = 2, pad = s.n))
            if k != s.sz_subspace
                println(io)
            end
        end
    else # never print more than 4 elements
        println(io, " "^(indent + 1), string(s[1]; base = 2, pad = s.n))
        println(io, " "^(indent + 1), string(s[2]; base = 2, pad = s.n))
        println(io, " "^(indent + 1), "⋮")
        println(io, " "^(indent + 1), string(s[end-1]; base = 2, pad = s.n))
        print(io, " "^(indent + 1), string(s[end]; base = 2, pad = s.n))
    end
end
