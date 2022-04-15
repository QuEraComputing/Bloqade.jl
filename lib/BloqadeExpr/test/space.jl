using Test
using BloqadeExpr

@testset "subspace" begin
    @test_throws ArgumentError space = Subspace(8, [1, 3, 5, 300])

    space = Subspace(8, [1, 3, 5, 10])
    @test space.map[3] == 2
    @test space.map[10] == 4

    state = rand(15)
    @test state[space] == state[[2, 4, 6, 11]]

    show(stdout, MIME"text/plain"(), space)
end
