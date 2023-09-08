using Test
using BloqadeKrylov
using LinearAlgebra
using SparseArrays
using BloqadeExpr



@testset "internal _mulp!" begin

    A = sprand(10, 10, 0.5)
    st = rand(10)

    G = A
    # testing p:
    for p in 1:5
        tar = similar(st)
        BloqadeKrylov._mulp!(tar, A, st, p)

        src = G*st

        @test tar ≈ src
        G *= A
    end

    
end


@testset "internal  _sign_roundup(X::T)" begin

    a = 3.34
    b = -3.24
    c = 0.0
    d = -0.0
    e = 1+1im
    f = -2+1im


    @test BloqadeKrylov._sign_roundup(a) == 1
    @test BloqadeKrylov._sign_roundup(b) == -1
    @test BloqadeKrylov._sign_roundup(c) == 1
    @test BloqadeKrylov._sign_roundup(d) == 1
    @test BloqadeKrylov._sign_roundup(e) == (1+1im)/√2
    @test BloqadeKrylov._sign_roundup(f) == (-2+1im)/√5
    
end


@testset "internal  _check_cols_parallel (REAL)" begin

    A = randn(10,10)
    A = BloqadeKrylov._sign_roundup.(A)

    @test BloqadeKrylov._check_cols_parallel(A,A) == true


    
end

@testset "internal  _check_cols_parallel (Complex)" begin

    A = randn(ComplexF64, 10,10)
    A = BloqadeKrylov._sign_roundup.(A)


    @test BloqadeKrylov._check_cols_parallel(A,A) == true
    
end



@testset "internal  _check_vecs_parallel (Complex)" begin

    A = randn(ComplexF64, 10,10)
    A = BloqadeKrylov._sign_roundup.(A)

    for i in 1:10
        for j in 1:10
            if i == j
                @test BloqadeKrylov._check_vecs_parallel(A[:,i],A[:,j]) == true
            end
        end
    end


    
end

@testset "onenormest impl" begin

    ## real
    A = [1. 0. 0. ; 5. 8. 2. ; 0 -1 0]
    
    nrm = LinearAlgebra.opnorm(A,1)
    nrmest,_,_ = BloqadeKrylov._onenormest_impl(A, transpose(A), 1)

    @test nrm ≈ nrmest


    A = [1. 0. 0. ; 5. 8. 2. ; 0 -1 0]
    
    nrm = LinearAlgebra.opnorm(A*A,1)
    nrmest,_,_ = BloqadeKrylov._onenormest_impl(A, transpose(A), 2)

    @test nrm ≈ nrmest

    ## complex
    A = [1+1im 2+0.5im 3.; 0. 6.7im -3.2-1im; 3.14 0.7+0.55im 0.29im+3.2]
    
    nrm = LinearAlgebra.opnorm(A,1)

    nrmest, _ , _ = BloqadeKrylov._onenormest_impl(A, adjoint(A), 1)

    @test nrm ≈ nrmest

    A = [1+1im 2+0.5im 3.; 0. 6.7im -3.2-1im; 3.14 0.7+0.55im 0.29im+3.2]
    
    nrm = LinearAlgebra.opnorm(A*A,1)

    nrmest, _ , _ = BloqadeKrylov._onenormest_impl(A, adjoint(A), 2)

    @test nrm ≈ nrmest


end


@testset "onenormest, API" begin

    ## real
    A = [1. 0. 0. ; 5. 8. 2. ; 0 -1 0]
    
    nrm = LinearAlgebra.opnorm(A,1)
    nrmest = BloqadeKrylov.onenormest(A, 1)

    @test nrm ≈ nrmest


    A = [1. 0. 0. ; 5. 8. 2. ; 0 -1 0]
    
    nrm = LinearAlgebra.opnorm(A*A,1)
    nrmest = BloqadeKrylov.onenormest(A, 2)

    @test nrm ≈ nrmest

    ## complex
    A = [1+1im 2+0.5im 3.; 0. 6.7im -3.2-1im; 3.14 0.7+0.55im 0.29im+3.2]
    
    nrm = LinearAlgebra.opnorm(A,1)

    nrmest = BloqadeKrylov.onenormest(A, 1)

    @test nrm ≈ nrmest

    A = [1+1im 2+0.5im 3.; 0. 6.7im -3.2-1im; 3.14 0.7+0.55im 0.29im+3.2]
    
    nrm = LinearAlgebra.opnorm(A*A,1)

    nrmest = BloqadeKrylov.onenormest(A, 2)

    @test nrm ≈ nrmest
    

end