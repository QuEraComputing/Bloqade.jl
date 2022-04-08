struct ComplexArray{T <: Real, N, AT <: AbstractArray{T}} <: AbstractArray{Complex{T}, N}
    storage::AT

    function ComplexArray{T, N}(::UndefInitializer, dims::Dims) where {T, N}
        storage = Array{T, length(dims)+1}(undef, (dims..., 2))
        new{T, N, typeof(storage)}(storage)
    end

    function ComplexArray{T, N}(storage::AbstractArray{T}) where {T, N}
        ndims(storage) == N+1 || throw(ArgumentError("storage dimension not match"))
        new{T, N, typeof(storage)}(storage)
    end
end

ComplexArray(storage::AbstractArray{T, N}) where {T, N} = ComplexArray{T, N-1}(storage)

Base.ndims(X::ComplexArray) = ndims(X.storage) - 1
Base.size(X::ComplexArray, idx::Int) = size(X.storage, idx)
Base.size(X::ComplexArray) = ntuple(idx->size(X, idx), ndims(X))

function Base.getindex(X::ComplexVector, idx::Int)
    @show "aaa"
    return X.storage[idx, 1] + im * X.storage[idx, 2]
end

function Base.getindex(X::ComplexArray, idx::Int...)
    return X.storage[idx..., 1] + im * X.storage[idx..., 2]
end

Adapt.@adapt_structure ComplexArray

function Adapt.adapt_storage(::Type{<:ComplexArray}, xs::AbstractArray{Complex{T}}) where T
    return ComplexArray(cat(real(xs), imag(xs); dims=ndims(xs)+1))
end
