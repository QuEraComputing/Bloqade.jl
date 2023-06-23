using Test
using BloqadeKrylov
using SparseArrays

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

#=
@testset "develop, please remove before pub" begin

    A = rand(10, 10)

    println( BloqadeKrylov.get_optimal_sm(1, A) )
    
end 
=#