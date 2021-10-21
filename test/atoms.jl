using Test
using LightGraphs
using RydbergEmulator
using Random

@testset "interfaces" begin
    atom = RydAtom([1, 2])
    @test atom.loc == (1, 2)
    @test atom[1] == 1
    @test atom[2] == 2
    atoms = rand_atoms(10, 0.2)
    @test all(ndims.(atoms) .== 2)

    atom = RydAtom((1, 2, 3))
    @test length(atom) == 3
    @test ndims(atom) == 3
    @test eltype(atom) === Int
end

@testset "canonical order" begin
    prev = -1
    for each in square_lattice(10, 0.8)
        @test each[1] >= prev
        prev = each[1]
    end
end

@test_logs (:warn, "graph has empty edges, creating a subspace contains the entire fullspace, consider using a full space register.") begin
    Random.seed!(123)
    blockade_subspace(SimpleGraph(5))
end
