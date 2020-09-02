using Test
using RydbergEmulator
using Random

atom = RydAtom([1, 2])
@test atom.loc == (1, 2)
@test atom[1] == 1
@test atom[2] == 2
atoms = rand_atoms(10, 0.2)
@test all(ndims.(atoms) .== 2)

@testset "canonical order" begin
    prev = -1
    for each in square_lattice(10, 0.8)
        @test each[1] >= prev
        prev = each[1]
    end    
end

@test_logs (:warn, "graph has empty edges, creating a subspace contains the entire fullspace, consider using a full space register.") begin
    Random.seed!(123)
    blockade_subspace(rand_atoms(4, 0.2))
end
