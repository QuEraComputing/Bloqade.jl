# bypass getindex of default broadcast
using Base.Broadcast: BroadcastStyle, ArrayStyle, AbstractArrayStyle, Broadcasted, DefaultArrayStyle

struct ComplexArrayStyle{S, N} <: AbstractArrayStyle{N} end

function ComplexArrayStyle{S, M}(::Val{N}) where {S, M, N}
    T = S <: AbstractArrayStyle{M} ? typeof(S(Val{N}())) : S
    return ComplexArrayStyle{T, N}()
end

function Broadcast.BroadcastStyle(::Type{AT}) where {T, N, ST, AT <: ComplexArray{T, N, ST}}
    return ComplexArrayStyle{typeof(BroadcastStyle(ST)), N}()
end

function Base.similar(bc::Broadcasted{ComplexArrayStyle{S, N}}, ::Type{Complex{T}}) where {S <: DefaultArrayStyle, N, T}
    return similar(ComplexArray{T}, axes(bc))
end

function Base.copyto!(dest::ComplexArray, bc::Broadcasted{<:ComplexArrayStyle})
    # fast path for operation only have `ComplexArray`
    # and real numbers
    if all_complex_array_type(bc.args...)
        args = unwarp_complex_array_type(bc.args...)
        broadcast!(bc.f, dest.storage, args...)
        return dest
    else
        return copyto!(dest, convert(Broadcasted{Nothing}, bc))
    end
end

# based on StructArrays, help aliasing analysis
Base.dataids(S::ComplexArray) = Base.dataids(S.storage)

unwarp_complex_array_type() = ()

function unwarp_complex_array_type(x, xs...)
    return (x, unwarp_complex_array_type(xs...)...)
end

function unwarp_complex_array_type(x::Complex, xs...)
    error("complex number broadcast is not supported for ComplexArray")
end

function unwarp_complex_array_type(x::AbstractArray, xs...)
    error("mixed array type broadcast is not supported for ComplexArray")
end

function unwarp_complex_array_type(x::ComplexArray, xs...)
    return (x.storage, unwarp_complex_array_type(xs...)...)
end

all_complex_array_type(x) = false
all_complex_array_type(x, xs...) = false
all_complex_array_type(::ComplexArray) = true
function all_complex_array_type(x::ComplexArray, xs...)
    all_complex_array_type(xs...)
end
