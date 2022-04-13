using YaoSubspaceArrayReg, Test

@testset "inplace reg operations" begin
    s = blockade_subspace([(0.0, 0.0), (0.6, 0.0), (1.2, 0.0)], 1.0)
    a = rand_state(s)
    regadd!(a, b)
end