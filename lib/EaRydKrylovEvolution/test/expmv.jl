using Test
using EaRydKrylovEvolution: expmv!

@testset "expmv" begin
    N = 400
    A = randn(ComplexF64, N, N)
    v = randn(ComplexF64, N)
    w = zeros(ComplexF64, N)
    A = A + A'
    expmv!(w, 10im, A, v)
    v3 = exp(10im*A) * v
    @test isapprox(w, v3; atol=1e-6)
end
