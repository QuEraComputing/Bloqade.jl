function Base.show(io::IO, ::MIME"text/plain", x::ComplexArray{T, N, <:CuArray}) where {T, N}
    summary(io, x)
    println(io)
    Base.print_array(io, adapt(Array, x))
end

# make CUDA work
function LinearAlgebra.norm(x::ComplexArray{T, N, <:CuArray}, p::Real=2) where {T, N}
    return @views norm(x.storage[:, 1] + x.storage[:, 2]*im, p)
end

# TODO: move this to YaoArrayRegister
function LinearAlgebra.normalize!(r::ArrayReg)
    return normalize!(statevec(r))
end

function Base.similar(A::ComplexArray{T, N, <:CuArray}, ::Type{S}, dims::Dims) where {T, N, S}
    return ComplexArray(similar(A.storage, real(S), (dims..., 2)))
end

function LinearAlgebra.rmul!(A::ComplexArray{T, N, <:CuArray{T}}, s::Complex) where {T, N}
    function kernel_fun(S, s)
        x, y = real(s), imag(s)
        im_stride = last(strides(S))
        idx = (blockIdx().x-1) * blockDim().x + threadIdx().x
        idx ≤ length(A) || return

        X_i = S[idx] * x - S[idx+im_stride] * y
        Y_i = S[idx] * y + S[idx+im_stride] * x
        S[idx          ] = X_i
        S[idx+im_stride] = Y_i
        return
    end

    total = length(A)
    kernel = @cuda launch=false kernel_fun(A.storage, s)
    threads = min(total, CUDA.maxthreads(kernel))
    blocks = cld(total, threads)
    kernel(A.storage, s; threads, blocks)
    return A
end
