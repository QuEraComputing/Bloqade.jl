using CUDA.CUSPARSE: CuSparseMatrixCSC,
    CuSparseMatrixCSR,
    AbstractCuSparseMatrix
using EaRydKrylov: expmv, expmv!
using SparseArrays
using LinearAlgebra
using BenchmarkTools
using EaRydCUDA
using CUDA

n = 18
A = sprand(1<<n, 1<<n, 1e-3)
x = rand(1<<n)
y = similar(x)
dA = cu(A)
dx = cu(x)
dy = similar(dx)

@benchmark mul!($y, $A, $x)
@benchmark mul!($dy, $dA, $dx)

@benchmark expmv!($y, 1.0, $A, $x)
@benchmark CUDA.@sync expmv!($dy, 1.0, $dA, $dx)

@profview expmv!(y, 1.0, A, x)
@profview expmv!(dy, 0.1, dA, dx)

@benchmark expmv!($y, 0.01, $A, $x)
@benchmark CUDA.@sync expmv!($dy, 0.01, $dA, $dx)

