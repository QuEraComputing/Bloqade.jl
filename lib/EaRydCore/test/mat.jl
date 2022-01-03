using Test
using Random, Yao
using EaRydCore
using Graphs
using SparseArrays

@testset "mat" begin
    ss = independent_set_subspace(smallgraph(:petersen))
    U = zeros(ComplexF64, 1024, 1024)
    U[ss, ss] = rand_unitary(length(ss))
    H = U * U'
    for g in [put(10, 2=>X), put(10, 3=>Rx(0.4)), put(10, 2:6=>matblock(rand_unitary(32)))*2, control(10, (3,-5), (2,7)=>matblock(rand_unitary(4))),
            chain(matblock(U), matblock(U)), igate(10), subroutine(10, kron(X,X), (5,2)), control(10, (3,-5), (2,7)=>matblock(rand_unitary(4)))',
            put(10, 2=>X) + matblock(rand_unitary(1024)), time_evolve(matblock(H), 0.3)
        ]
        M = SparseMatrixCSC(mat(g))[ss, ss]
        @test M ≈ mat(g, ss)
    end
end