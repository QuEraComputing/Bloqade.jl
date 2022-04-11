function LinearAlgebra.mul!(y::ComplexVector, A::AbstractVecOrMat{<:Real}, x::ComplexVector, alpha::Number, beta::Number)
    mul!(y.storage, A, x.storage, alpha, beta)
    return y
end

function LinearAlgebra.mul!(y::ComplexVector, A::Diagonal{<:Real}, x::ComplexVector, alpha::Number, beta::Number)
    mul!(y.storage, A, x.storage, alpha, beta)
    return y
end

# (X - im*Y)(A + im*B)
# sum_i (X_i - im Y_i)(A_i + im*B_i)
# sum_i X_i A_i + X_i * B_i * im - Y_i * A_i * im + Y_i * B_i
# (X⋅A+Y⋅B) + (X⋅B-Y⋅A)im
function LinearAlgebra.dot(x::ComplexVector, y::ComplexVector)
    @views begin
        X = x.storage[:, 1]
        Y = x.storage[:, 2]
        A = y.storage[:, 1]
        B = y.storage[:, 2]            
    end
    return Complex(dot(X, A) + dot(Y, B), dot(X, B)-dot(Y, A))
end

# (A+im*B)*a + (X+im*Y)
# A*a + im*B*a + X + im*Y
# X <- (A*a + X)
# Y <- (B*a + Y)
function LinearAlgebra.axpy!(a::Real, x::ComplexVector, y::ComplexVector)
    @views begin
        A = x.storage[:, 1]
        B = x.storage[:, 2]
        X = y.storage[:, 1]
        Y = y.storage[:, 2]
    end
    axpy!(a, A, X)
    axpy!(a, B, Y)
    return y
end

# (A+im*B)*(a+b*im) + (X+im*Y)
# A*a + A*b*im + im*B*a - B*b + X + im*Y
# X <- (A*a - B*b + X)
# Y <- (A*b + B*a + Y)
#
# X <- A * a + X
# X <- -b * B + (A * a + X)
#
# Y <- A*b + Y
# Y <- B*a + (A*b + Y)
function LinearAlgebra.axpy!(alpha::Complex, x::ComplexVector, y::ComplexVector)
    @views begin
        A = x.storage[:, 1]
        B = x.storage[:, 2]
        X = y.storage[:, 1]
        Y = y.storage[:, 2]
    end
    a, b = real(alpha), imag(alpha)
    axpy!(a, A, X)
    axpy!(-b, B, X)
    axpy!(b, A, Y)
    axpy!(a, B, Y)
    return y
end

function Base.reshape(X::ComplexArray, dims::Dims)
    return ComplexArray(reshape(X.storage, (dims..., 2)))
end

function LinearAlgebra.rmul!(X::ComplexArray, s::Real)
    rmul!(X.storage, s)
    return X
end

# (X + im * Y) * (x + im * y)
# X * x + im * X * y + im * Y * x - Y * y
# X * x - Y * y + im * X * y + im * Y * x
# X <- X * x - Y * y
# Y <- X * y + Y * x
function LinearAlgebra.rmul!(A::ComplexArray{T, N}, s::Complex) where {T, N}
    x, y = real(s), imag(s)
    im_stride = last(strides(A.storage))
    @inbounds for i in 1:length(A)
        X_i = A.storage[i] * x - A.storage[i+im_stride] * y
        Y_i = A.storage[i] * y + A.storage[i+im_stride] * x
        A.storage[i          ] = X_i
        A.storage[i+im_stride] = Y_i
    end
    return A
end
