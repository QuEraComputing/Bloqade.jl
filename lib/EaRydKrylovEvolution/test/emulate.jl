using Test
using EaRydKrylovEvolution
using EaRydKrylovEvolution: emulate_routine!, nsites, get_space, PrecisionAdaptor, adapt
using SparseArrays
using LinearAlgebra

@testset "emulate_routine" begin
    atoms = square_lattice(10, 0.8)
    t = 0.1
    h = rydberg_h(atoms; Ω=0.1, Δ=0.2)

    @testset "fullspace" begin
        cache = SparseMatrixCSC(h)
        result = statevec(emulate_routine!(zero_state(10), t, h, cache))
        answer = exp(-im * t * Matrix(cache)) * statevec(zero_state(10))
        @test result ≈ answer
    end

    @testset "subspace" begin
        space = blockade_subspace(atoms)
        cache = SparseMatrixCSC(h, space)
        result = statevec(emulate_routine!(zero_state(space), t, h, cache))
        answer = exp(-im * t * Matrix(cache)) * statevec(zero_state(space))
        @test result ≈ answer
    end

    @testset "real layout" begin
        space = blockade_subspace(atoms)
        cache = SparseMatrixCSC(h, space)
        @test_throws ErrorException emulate_routine!(zero_state(space, RealLayout()), t, h, cache)
        # result = statevec(emulate_routine!(zero_state(space, RealLayout()), t, h, cache))
        # answer = exp(-im * t * Matrix(cache)) * statevec(zero_state(space, RealLayout()))
        # @test result ≈ answer
    end

    for space in [FullSpace(), blockade_subspace(atoms)]
        cache = KrylovEmulationCache{Float64, Cint}(h, FullSpace())
        @test nnz(cache.H) == nnz(SparseMatrixCSC(h))
        @test typeof(cache.H) === SparseMatrixCSC{Float64, Cint}
    end
end

function naive_discrete_evolve(P, reg, ts, hs)
    st = statevec(adapt(PrecisionAdaptor(P), reg))
    for (t, h) in zip(ts, hs)
        H = Matrix(SparseMatrixCSC(h, get_space(reg)))
        st = exp(-im * t * H) * st
    end
    return st
end

@testset "KrylovEvolution" begin
    atoms = square_lattice(10, 0.8)
    space = blockade_subspace(atoms)

    @testset "QAOA P=$P reg=$(nameof(typeof(reg)))" for reg in [zero_state(10), zero_state(space)],
            P in [Float32, Float64]

        durations = rand(5)
        hs = [rydberg_h(atoms; Ω, Δ) for (Ω, Δ) in zip(rand(5), rand(5))]

        evolve = KrylovEvolution{P}(copy(reg), durations, hs)
        @test eltype(evolve.durations) === P
        @test eltype(evolve.cache.H) === P
        @test eltype(evolve.reg.state) === Complex{P}

        emulate!(evolve)
        @test statevec(evolve.reg) ≈ naive_discrete_evolve(P, reg, durations, hs) rtol=sqrt(eps(P))
    end

    @testset "discretize continuous" for total_time in [0.1, 0.1f0]
        P = typeof(total_time)
        h = rydberg_h(atoms; Ω=sin, Δ=cos)
        durations, hs = trotterize(total_time, h, nsteps=10)
        evolve = KrylovEvolution(durations, hs)

        @test eltype(durations) == P
        @test evolve.reg isa ArrayReg
        @test evolve.durations == durations
        @test eltype(evolve.reg.state) === Complex{eltype(durations)}

        emulate!(evolve)
        reg = zero_state(nsites(first(hs)))
        target = naive_discrete_evolve(P, reg, durations, hs)
        @test statevec(evolve.reg) ≈ target rtol=sqrt(eps(P))
    end

    @testset "measure observables" for total_time in [0.1, 0.1f0]
        total_time = 0.1
        P = typeof(total_time)
        h = rydberg_h(atoms; Ω=sin, Δ=cos)
        durations, hs = trotterize(total_time, h, nsteps=10)
        evolve = KrylovEvolution(durations, hs)
        observables = Float64[]
        for (_, reg, _, _) in evolve
            push!(observables, expect(put(10, 1=>X), reg))
        end

        st = statevec(zero_state(nsites(h)))
        target = Float64[]
        for (t, h) in zip(durations, hs)
            H = Matrix(SparseMatrixCSC(h))
            st = exp(-im * t * H) * st
            push!(target, expect(put(10, 1=>X), ArrayReg(st)))
        end

        @test observables ≈ target
    end
end
