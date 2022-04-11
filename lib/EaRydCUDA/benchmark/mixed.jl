# mixed normal/complex type SpMV benchmark
using CUDA
using CUDA.CUSPARSE
using LinearAlgebra
using SparseArrays
using BenchmarkTools

function benchmark_suite(::Type{T}, n::Int) where T
    S = sprand(T, 1<<n, 1<<n, 1e-3)
    X = rand(Complex{T}, 1<<n)
    Y = similar(X)

    dS = CuSparseMatrixCSC(S)
    dX = CuVector(X)
    dY = CuVector(Y)

    cpu = @benchmark mul!($Y, $S, $X)
    cuda = @benchmark CUDA.@sync mul!($dY, $dS, $dX)

    A = sprand(Complex{T}, 1<<n, 1<<n, 1e-3)
    dA = CuSparseMatrixCSC(A)
    cuda_complex = @benchmark CUDA.@sync mul!($dY, $dA, $dX)

    minimum(cpu).time, minimum(cuda).time, minimum(cuda_complex).time
end

report = (cpu=Float64[], cuda=Float64[], cuda_complex=Float64[])
for n in 10:15
    @info "benchmarking..." n
    cpu, cuda, cuda_complex = benchmark_suite(Float32, n)
    push!(report.cpu, cpu)
    push!(report.cuda, cuda)
    push!(report.cuda_complex, cuda_complex)
end

report.cpu ./ report.cuda
report.cpu ./ report.cuda_complex
