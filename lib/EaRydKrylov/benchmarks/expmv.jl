# benchmark expmv performance for linear map vs sparse matrix
# conclusion: they have similar performance, lienar map is slightly
# slower

using Test
using EaRydKrylov
using EaRydExpr
using EaRydExpr: Hamiltonian
using BenchmarkTools

function benchmark(n)
    atom = [(i, ) for i in 1:n]
    h = rydberg_h(atom, Ω=0.5, C=109.23)
    H = Hamiltonian(Float64, h)

    M = sum(zip(H.fs, H.ts)) do (f, t)
        f(0.1) * t
    end

    st = rand(ComplexF64, 1<<n)
    w = zeros(ComplexF64, 1<<n)

    new = @benchmark EaRydKrylov.expmv!($w, 0.1im, $(H(0.1)), $st)
    old = @benchmark EaRydKrylov.expmv!($w, 0.1im, $(M), $st)
    return minimum(new).time, minimum(old).time
end

results = map(benchmark, [5, 8, 10, 15])
