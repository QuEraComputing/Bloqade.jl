using Test
using BloqadeExpr
using BloqadeKrylov


@testset "ValHamiltonian" begin

    Ham = BloqadeExpr.Hamiltonian(Float64, SumOfX(1, sin) + SumOfZ(1,cos))


    t = 0.645
    StepHam = Ham(0.645)
    ValHam = BloqadeKrylov.Val(StepHam)


    # check coefficents:
    for (i,f) in enumerate(StepHam.h.fs)
        @test f(t) == ValHam.fvals[i]
    end

    @test size(StepHam) == size(ValHam)

    @test StepHam.h === ValHam.h
    @test ValHam.h === Ham


    # check basic algos :+
    AddVHam = ValHam + ValHam
    @test AddVHam.h === ValHam.h
    @test AddVHam.h === Ham

    for (i,f) in enumerate(StepHam.h.fs)
        @test AddVHam.fvals[i] == f(t) + f(t)
    end

    # check basic algos :-
    SubVHam = ValHam - ValHam
    @test SubVHam.h === ValHam.h
    @test SubVHam.h === Ham

    for (i,f) in enumerate(StepHam.h.fs)
        @test SubVHam.fvals[i] == f(t) - f(t)
    end

    # check basic algos :*
    MulVHam = 0.5*ValHam 
    @test MulVHam.h === ValHam.h
    @test MulVHam.h === Ham

    for (i,f) in enumerate(StepHam.h.fs)
        @test MulVHam.fvals[i] == 0.5*f(t)
    end


end