using CUDA.CUSPARSE: CuSparseMatrixCSC,
    CuSparseMatrixCSR,
    AbstractCuSparseMatrix
using EaRydKrylov: expmv, expmv!
using SparseArrays
using LinearAlgebra
using BenchmarkTools
using ComplexArrays
using EaRydCUDA
using CUDA
CUDA.allowscalar(false)

T = ComplexF32
n = 15
x = rand(T, 1<<n)
y = rand(T, 1<<n)
A = sprand(T, 1<<n, 1<<n, 1e-3)
cx = adapt(ComplexArray, x)
cy = adapt(ComplexArray, y)
dx = adapt(CuArray, x)
dy = adapt(CuArray, y)
dA = adapt(CuArray, A)
dcx = adapt(CuArray, cx)
dcy = adapt(CuArray, cy)

# ComplexVector is 4x slower
@benchmark dot($x, $y)
@benchmark CUDA.@sync dot($dx, $dy)
@benchmark dot($cx, $cy)
@benchmark CUDA.@sync dot($dcx, $dcy)


@benchmark axpy!(1.0, $x, $y)
@benchmark CUDA.@sync axpy!(1.0, $dx, $dy)
@benchmark axpy!(1.0, $cx, $cy)
@benchmark CUDA.@sync axpy!(1.0, $dcx, $dcy)

A = sprand(real(eltype(x)), 1<<n, 1<<n, 1e-3)
dA = adapt(CuArray, A)

mul_cpu_normal = @benchmark mul!($y, $A, $x)
@benchmark mul!($dy, $dA, $dx) # doesn't work

mul_cpu_complex = @benchmark mul!($cy, $A, $cx) # ~2x slower
mul_cuda_complex = @benchmark CUDA.@sync mul!($dcy, $dA, $dcx)

minimum(mul_cpu_normal).time/minimum(mul_cpu_complex).time
minimum(mul_cpu_normal).time/minimum(mul_cuda_complex).time