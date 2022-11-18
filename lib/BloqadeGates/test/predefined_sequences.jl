using BloqadeGates
using BloqadeExpr
using BloqadeODE, BloqadeKrylov
using Yao
using Test

using BloqadeGates: two_level_indices

function test_two_sites_pulse_sequence(reg, sequence; op = SWAP, atol = 1e-6)
    goal = copy(state(reg))
    M = nothing
    for p in sequence
        mh = mat(p)
        isnothing(M) ? (M = mh) : (M = mh*M)
        goal = mh * goal
        reg |> p
        fdlt = abs((goal' * state(reg))[1,1])
        @test isapprox(fdlt, 1; atol)
    end
    if nqudits(reg) == 2
        ids = two_level_indices(2)
        @test operator_fidelity(matblock(M[ids, ids]), op) > 1 - atol
    end
end

@testset "Local pulses" begin
    @testset "C_{$n}Z" for n = 1:5
        atoms = [(cos(2π*i/n), sin(2π*i/n)) for i = 1:n]
        if n == 1
            @test_throws AssertionError local_CkZ(atoms, [i for i = 2:n], [1])
        else
            seq = local_CkZ(atoms, [i for i = 2:n], [1])
            m = mat(seq)
            ids = two_level_indices(n)
            m = m[ids, ids]
            @test operator_fidelity(matblock(m), control(n, collect(2:n), 1=>Z)) ≈ 1
        end
    end
    
    @testset "C_{$n}NOT" for n = 1:5
        atoms = [(cos(2π*i/n), sin(2π*i/n)) for i = 1:n]
        if n == 1
            @test_throws AssertionError local_CkZ(atoms, [i for i = 2:n], [1])
        else
            seq = local_CkNOT(atoms, [i for i = 2:n], [1])
            m = mat(seq)
            ids = two_level_indices(n)
            m = m[ids, ids]
            @test operator_fidelity(matblock(m), control(n, collect(2:n), 1=>X)) ≈ 1
        end
    end
end
