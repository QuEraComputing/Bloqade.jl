using Test
using BloqadeKrylov
using SparseArrays
using BloqadeLattices
using BloqadeWaveforms
using BloqadeExpr: Hamiltonian 
using BloqadeExpr
using Yao

@testset "expm_multiply_impl Real, dense A" begin

    A = rand(10, 10)
    st = rand(10)


    vold = BloqadeKrylov.expmv(1, A, st)


    vnew = similar(st)
    st_before = deepcopy(st)
    BloqadeKrylov._expm_multiply_impl!(vnew, 1, A, st, 0, 1, 40, 2^-53)

    # checking st is not modefied
    @test st_before == st
    
    #  checking results against expmv
    @test vnew ≈ vold
end

@testset "expm_multiply_impl complex, sparse A" begin

    A = sprand(ComplexF64,10, 10, 0.6)
    st = rand(ComplexF64,10)


    vold = BloqadeKrylov.expmv(1, A, st)


    vnew = similar(st)
    st_before = deepcopy(st)
    BloqadeKrylov._expm_multiply_impl!(vnew, 1, A, st, 0, 1, 40, 2^-53)

    # checking st is not modefied
    @test st_before == st
    
    #  checking results against expmv
    @test vnew ≈ vold

end


@testset "expm_multiply Real, dense A" begin

    A = rand(10, 10)
    st = rand(10)


    vold = BloqadeKrylov.expmv(2.55, A, st)


    vnew = similar(st)
    st_before = deepcopy(st)
    BloqadeKrylov.expm_multiply!(vnew, 2.55, A, st)

    # checking st is not modefied
    @test st_before == st
    
    #  checking results against expmv
    @test vnew ≈ vold

end

@testset "expm_multiply complex, sparse A" begin

    A = sprand(ComplexF64,10, 10, 0.6)
    st = rand(ComplexF64,10)


    vold = BloqadeKrylov.expmv(2.55, A, st)


    vnew = similar(st)
    st_before = deepcopy(st)
    BloqadeKrylov.expm_multiply!(vnew, 2.55, A, st)

    # checking st is not modefied
    @test st_before == st
    
    #  checking results against expmv
    @test vnew ≈ vold

end

@testset "expm_multiply, Rydberg" begin

    atoms = generate_sites(ChainLattice(), 4, scale = 6.1)
    clocks = [0.0, 0.5, 0.8, 1.1, 1.5, 1.6]
    wf = piecewise_constant(clocks = clocks, values = [0.0, 2.1, 2.1, 1.5, 0.0])
    h = rydberg_h(atoms; Ω = wf)

    Ham = Hamiltonian(Float64, h)
    #print(Ham)
    Ht = Ham(0.6)
    reg = zero_state(length(atoms))
    state = statevec(reg)


    vold = BloqadeKrylov.expmv(0.05, to_matrix(Ht), state)

    println("main")
    println(vold)
    vs = similar(state)
    BloqadeKrylov.expm_multiply!(vs, 0.05, to_matrix(Ht), state)
    #println(vs)
    #vnew = BloqadeKrylov.expm_multiply(2.55, to_matrix(Ht), state)

    #println(vnew)
    @test vold ≈ vs


end


#=
@testset "develop, please remove before pub" begin

    A = rand(10, 10)

    println( BloqadeKrylov.get_optimal_sm(1, A) )
    
end 
=#