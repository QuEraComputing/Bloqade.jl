using BloqadeLattices
using Test

@testset "rydberg_interaction_matrix" begin
    pbc = parallelepiped_region(ChainLattice(),(5,);pbc=true)
    obc =  parallelepiped_region(ChainLattice(),(5,);pbc=false)

    V_pbc = zeros(5,5)
    V_pbc[1,2] = V_pbc[2,3] = V_pbc[3,4] = V_pbc[4,5] = 1
    V_pbc[1,3] = V_pbc[2,4] = V_pbc[3,5] = 1/2^6
    V_pbc[1,4] = V_pbc[2,5] = 1/2^6
    V_pbc[1,5] = 1
    
    @test V_pbc == rydberg_interaction_matrix(pbc,1)

    V_obc = zeros(5,5)
    V_obc[1,2] = V_obc[2,3] = V_obc[3,4] = V_obc[4,5] = 1
    V_obc[1,3] = V_obc[2,4] = V_obc[3,5] = 1/2^6
    V_obc[1,4] = V_obc[2,5] = 1/3^6
    V_obc[1,5] = 1/4^6

    @test V_obc == rydberg_interaction_matrix(obc,1)

end