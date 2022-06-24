using Test
using BloqadeExpr

@testset "subspace" begin
    @test_throws ArgumentError space = Subspace(8, [1, 3, 5, 300])

    # 1, 3, 5, 10 are full space indexes
    # map[fullspace_index] = subspace index
    space = Subspace(8, [1, 3, 5, 10])
    @test space.map[3] == 2
    @test space.map[10] == 4

    for pair in space
        @test pair.first in keys(space)
        @test pair.second in values(space)
    end

    @test space[10] == 4

    @test space[bit"1010"] == 4

    @test keys(space) == Set([1, 3, 5, 10])

    @test sort(collect(values(space))) == [1, 2, 3, 4]

    @test haskey(space, 1)
    
    copied_space = copy(space)
    @test copied_space.nqubits == space.nqubits
    @test copied_space.map == space.map
    @test copied_space.subspace_v == space.subspace_v

    @test vec(space) == [1, 3, 5, 10]

    state = rand(15)
    @test state[space] == state[[2, 4, 6, 11]]
    
    # non-compact print
    show(stdout, MIME"text/plain"(), space)

    print("\n\n")

    # compact print
    space = Subspace(6, [i for i in 1:40])
    show(IOContext(stdout, :limit => true), MIME"text/plain"(), space)

    # all contents equal
    space1 = Subspace(6, [i for i in 1:40])
    space2 = Subspace(6, [i for i in 1:40])
    @test space1 == space2

    # nqubits mismatch
    space2 = Subspace(7, [i for i in 1:40])
    @test !(space1 == space2)

    # subspace values length mismatch
    space2 = Subspace(6, [i for i in 1:20])
    @test !(space1 == space2)

    # map mismatch
    space2 = Subspace(6, [i for i in 5:44])
    @test !(space1 == space2)

end
