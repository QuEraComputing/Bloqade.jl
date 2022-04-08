struct ComplexArray{T <: Real, N, AT <: AbstractArray{T}} <: AbstractArray{Complex{T}, N}
    storage::AT

    function ComplexArray{T, N}(::UndefInitializer, dims::Dims{N}) where {T, N}
        storage = Array{T, length(dims)+1}(undef, (dims..., 2))
        new{T, N, typeof(storage)}(storage)
    end

    function ComplexArray{T, N}(storage::AbstractArray{T}) where {T, N}
        ndims(storage) == N+1 || throw(ArgumentError("storage dimension not match"))
        new{T, N, typeof(storage)}(storage)
    end
end

const ComplexVector{T} = ComplexArray{T, 1}
const ComplexMatrix{T} = ComplexArray{T, 2}

ComplexArray{T}(::UndefInitializer, dims::Dims{N}) where {T, N} = ComplexArray{T, N}(undef, dims)
ComplexArray{T, N}(::UndefInitializer, dims::Int...) where {T, N} = ComplexArray{T, N}(undef, dims)
ComplexArray{T}(::UndefInitializer, dims::Int...) where {T, N} = ComplexArray{T}(undef, dims)
ComplexArray(storage::AbstractArray{T, N}) where {T, N} = ComplexArray{T, N-1}(storage)

Base.ndims(X::ComplexArray) = ndims(X.storage) - 1
Base.size(X::ComplexArray, idx::Int) = size(X.storage, idx)
Base.size(X::ComplexArray) = ntuple(idx->size(X, idx), ndims(X))

# NOTE: this is mainly for printing since Julia
# treat Vector printing same as Matrix printing
function Base.getindex(X::ComplexVector, idx::Int, others::Int...)
    all(isequal(1), others) || throw(BoundsError(X, (idx, others...)))
    return X.storage[idx, 1] + im * X.storage[idx, 2]
end

function Base.getindex(X::ComplexArray{T, N}, idx::Int) where {T, N}
    return Complex(X.storage[idx], X.storage[idx+last(strides(X.storage))])
end

function Base.getindex(X::ComplexArray{T, N}, idx::Vararg{Int, N}) where {T, N}
    return Complex(X.storage[idx..., 1], X.storage[idx..., 2])
end

function Base.setindex!(X::ComplexArray{T, N}, value::Complex{T}, idx::Int) where {T, N}
    X.storage[idx] = real(value)
    X.storage[idx+last(strides(X.storage))] = imag(value)
    return X
end

function Base.setindex!(X::ComplexArray{T, N}, value::Complex{T}, idx::Int, others::Int...) where {T, N}
    X.storage[idx, others..., 1] = real(value)
    X.storage[idx, others..., 2] = imag(value)
    return X
end

Base.IndexStyle(::Type{<:ComplexArray}) = IndexLinear()
function Base.similar(A::ComplexArray, ::Type{S}, dims::Dims) where {T, S <: Complex{T}}
    return ComplexArray{T}(undef, dims)
end

Adapt.@adapt_structure ComplexArray

function Adapt.adapt_storage(::Type{<:ComplexArray}, xs::AbstractArray{Complex{T}}) where T
    return copyto!(similar(ComplexArray{T}, size(xs)), xs)
end
