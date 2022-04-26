using BloqadeSchema
using Test
using Configurations

@testset "to_lattice with 1d chains" begin
    @test BloqadeSchema.to_lattice([1, 2, 3, 4, 5]) == BloqadeSchema.Lattice(;
        sites=[Tuple([1, 0]), Tuple([2, 0]), Tuple([3, 0]), Tuple([4, 0]), Tuple([5, 0])],
        filling=[1, 1, 1, 1, 1]
    )

    @test BloqadeSchema.to_lattice([]) == BloqadeSchema.Lattice(;
        sites=[],
        filling=[]
    )
end