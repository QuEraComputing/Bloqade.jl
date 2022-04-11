using CUDA.CUSPARSE: CuSparseMatrixCSC,
    CuSparseMatrixCSR,
    AbstractCuSparseMatrix
using EaRydKrylov: expmv, expmv!
using SparseArrays
using LinearAlgebra
using BenchmarkTools
using EaRydCUDA
using CUDA

CUDA.allowscalar(false)

# n = 15
# A = sprand(ComplexF32, 1<<n, 1<<n, 1e-3)
# x = rand(ComplexF32, 1<<n)
# y = similar(x)
# dA = cu(A)
# dx = cu(x)
# dy = similar(dx)

# @benchmark mul!($y, $A, $x)
# @benchmark mul!($dy, $dA, $dx)

# @benchmark expmv!($y, 1.0, $A, $x)
# @benchmark CUDA.@sync expmv!($dy, 1.0, $dA, $dx)

# @profview expmv!(y, 1.0, A, x)
# @profview expmv!(dy, 0.1, dA, dx)

# @benchmark expmv!($y, 0.01, $A, $x)
# @benchmark CUDA.@sync expmv!($dy, 0.01, $dA, $dx)

function benchmark_report(::Type{T}, range, t::Real=0.1) where T
    report = (
        blas_cpu=Float64[],
        blas_cuda=Float64[],
        expmv_cpu=Float64[],
        expmv_cuda=Float64[],
    )

    for n in range
        @info "benchmarking..." n
        blas_cpu, blas_cuda, expmv_cpu, expmv_cuda = benchmark_suite(T, n, t)
        push!(report.blas_cpu, blas_cpu)
        push!(report.blas_cuda, blas_cuda)
        push!(report.expmv_cpu, expmv_cpu)
        push!(report.expmv_cuda, expmv_cuda)
    end
    return report
end

function benchmark_suite(::Type{T}, n::Int, t::Real=0.1) where T
    A = sprand(T, 1<<n, 1<<n, 1e-3)
    x = rand(T, 1<<n)
    y = similar(x)
    dA = adapt(CuArray, A)
    dx = adapt(CuArray, x)
    dy = similar(dx)

    @info "arguments" y=typeof(y) A=typeof(A) x=typeof(x)
    blas_cpu = @benchmark mul!($y, $A, $x)
    @info "arguments" dy=typeof(dy) dA=typeof(dA) dx=typeof(dx)
    blas_cuda = @benchmark mul!($dy, $dA, $dx)

    expmv_cpu = @benchmark expmv!($y, $t, $A, $x)
    expmv_cuda = @benchmark CUDA.@sync expmv!($dy, $t, $dA, $dx)

    return minimum(blas_cpu).time,
        minimum(blas_cuda).time,
        minimum(expmv_cpu).time,
        minimum(expmv_cuda).time
end

report = Dict()
for T in [Float32, Float64, ComplexF32, ComplexF64]
    @info "benchmarking..." T
    report[T] = benchmark_report(T, 5:18, 0.1)
end

function generate_report(::Type{T}, report) where T
    report = report[T]
    blas_speedup = report.blas_cpu ./ report.blas_cuda
    expmv_speedup = report.expmv_cpu ./ report.expmv_cuda
    println("BLAS speedup $T:")
    display(blas_speedup)
    println()
    println("expmv speedup $T:")
    display(expmv_speedup)
end

for T in [Float32, Float64, ComplexF32, ComplexF64]
    generate_report(T, report)
end
