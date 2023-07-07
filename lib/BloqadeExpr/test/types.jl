using Test
using BloqadeExpr, YaoBlocks, YaoAPI, BitBasis
using SparseArrays
using LinearAlgebra


@testset "getindex" begin
    pb = rydberg_h([(5*randn(), 5*randn()) for i=1:5]; Ω=0.3, Δ=0.5)
    mpb = mat(pb)
    allpass = true
    for i=basis(pb), j=basis(pb)
        allpass &= pb[i, j] ≈ mpb[Int(i)+1, Int(j)+1]
    end
    @test allpass

    allpass = true
    for j=basis(pb)
        allpass &= vec(pb[:, j]) ≈ mpb[:, Int(j)+1]
        allpass &= vec(pb[j,:]) ≈ mpb[Int(j)+1,:]
        allpass &= vec(pb[:, EntryTable([j], [1.0+0im])]) ≈ mpb[:, Int(j)+1]
        allpass &= vec(pb[EntryTable([j], [1.0+0im]),:]) ≈ mpb[Int(j)+1,:]
    end
    @test allpass
end

@testset "isreal" begin
    @test isreal(SumOfN(nsites=2))
    @test isreal(SumOfZ(nsites=2))

    scale = 2 * SumOfX(nsites=3)
    @test isreal(scale)
    scale = im * SumOfX(nsites=3)
    @test !isreal(scale)

    # 3-level
    atoms = [(rand(), rand()) for i=1:3]
    @test !isreal(rydberg_h(atoms; Ω = 1.0, ϕ = sin, Δ = 1.0))
    @test !isreal(rydberg_h(atoms; ϕ = sin, Δ = 1.0))
    @test isreal(rydberg_h(atoms; Ω = sin, Δ = 1.0))

    @test !isreal(rydberg_h_3(atoms; Ω_r = 1.0, ϕ_r = sin, Ω_hf = 1.0, ϕ_hf = sin))
    @test !isreal(rydberg_h_3(atoms; Ω_r = 1.0, Ω_hf = 1.0, ϕ_hf = sin))
    @test isreal(rydberg_h_3(atoms; Ω_r = 1.0, Ω_hf = sin))
end

@testset "Hamiltonian Term Default Values" begin
    pauli_zs_term = SumOfZ(2)
    @test pauli_zs_term.Δ == 1
    
    number_ops_term = SumOfN(2)
    @test number_ops_term.Δ == 1

    pauli_xs_term = SumOfX(2)
    @test pauli_xs_term.Ω == 1
end

@testset "Qudit Numbers" begin
    @test YaoAPI.nqudits(PdPhase(1.0)) == 1
    @test YaoAPI.nqudits(PuPhase(1.0)) == 1
    @test YaoAPI.nqudits(SumOfZ(10)) == 10

    # 3-level
    @test YaoAPI.nqudits(PdPhase_01(π)) == 1
    @test YaoAPI.nqudits(PdPhase_1r(0.0)) == 1
    @test YaoAPI.nqudits(PuPhase_01(π)) == 1
    @test YaoAPI.nqudits(PuPhase_1r(0.0)) == 1
    @test YaoAPI.nqudits(SumOfZ_01(10)) == 10
end

@testset "Equality Operator" begin
    @test SumOfZ(10, 1.0) == SumOfZ(10, 1.0)
    @test !(SumOfZ(1, 1.0) == SumOfZ(2, 1.0))
    @test !(SumOfZ(5, 0.5) == SumOfZ(5, 0.3))

    h = rydberg_h([(5*randn(), 5*randn()) for i=1:5]; Ω=0.3, Δ=0.5)
    @test h == copy(h)
    h3 = rydberg_h_3([(5*randn(), 5*randn()) for i=1:5]; Ω_r=0.3, Δ_r=0.5, Ω_hf=0.3, Δ_hf=0.5)
    @test h3 == copy(h3)
end

@testset "Hamiltonian Size" begin
    hamiltonian = BloqadeExpr.Hamiltonian(Float64, SumOfX(6, sin) + SumOfZ(6, cos))
    @test size(hamiltonian) == (64, 64)
    @test size(hamiltonian, 2) == 64

    #if REPL is 32bit or 64bit Int for sparse will have differet size:
    if sizeof(Int) == 8
        @test storage_size(hamiltonian) == 6672
    else
        @test storage_size(hamiltonian) == 4872
    end
    step_hamiltonian = hamiltonian(0.1)
    @test size(step_hamiltonian) == (64, 64)
    @test size(step_hamiltonian, 2) == 64
end

@testset "Hamiltonian/SumOfLinop type" begin
    hamiltonian = BloqadeExpr.Hamiltonian(Float64, SumOfX(6, sin) + SumOfZ(6, cos))
    @test precision_type(hamiltonian) == Float64
    @test highest_type(hamiltonian) == Float64
    
    step_hamiltonian = hamiltonian(0.1)
    @test precision_type(step_hamiltonian) == Float64
    @test highest_type(step_hamiltonian) == Float64
end

@testset "SumOfLinop Norm" begin 
    hamiltonian = BloqadeExpr.Hamiltonian(Float64, SumOfZ(1, sin))
    step_hamiltonian = hamiltonian(π/2)
    @test LinearAlgebra.opnorm(step_hamiltonian,1) == 1.0
end


@testset "SumOfLinop " begin


    Ham = BloqadeExpr.Hamiltonian(Float64, SumOfX(1, sin) + SumOfZ(1,cos))


    t = 0.645
    StepHam = Ham(0.645)

<<<<<<< HEAD
    @test ishermitian(StepHam) == true
    

    # check coefficents:
    for (i,f) in enumerate(Ham.fs)
        @test f(t) == StepHam.fvals[i]
    end
    @test StepHam.ts === Ham.ts


    #=
    # check basic algos :+
    AddOp = StepHam + StepHam
    @test AddOp.ts === StepHam.ts
    @test AddOp.ts === Ham.ts

    for (i,f) in enumerate(Ham.fs)
        @test AddOp.fvals[i] == f(t) + f(t)
    end

    @test LinearAlgebra.is_hermitian(AddOp) == true

    # check basic algos :-
    SubVHam = StepHam - StepHam
    @test SubVHam.ts === StepHam.ts
    @test SubVHam.ts === Ham.ts

    for (i,f) in enumerate(Ham.fs)
        @test SubVHam.fvals[i] == f(t) - f(t)
    end

    @test LinearAlgebra.is_hermitian(SubVHam) == true
    =#

    # check basic algos :*
    MulVHam = 0.5*StepHam 
    @test MulVHam.ts === StepHam.ts
    @test MulVHam.ts === Ham.ts

    for (i,f) in enumerate(Ham.fs)
        @test MulVHam.fvals[i] == 0.5*f(t)
    end

    @test ishermitian(MulVHam) == true
    @test isskewhermitian(MulVHam) == false


     # check basic algos :*
     MulVHam2 = 0.5im*StepHam 
     @test MulVHam2.ts === StepHam.ts
     @test MulVHam2.ts === Ham.ts
 
     for (i,f) in enumerate(Ham.fs)
         @test MulVHam2.fvals[i] == 0.5im*f(t)
     end   

     @test ishermitian(MulVHam2) == false
     @test isskewhermitian(MulVHam2) == true


     # check basic algos :*
     MulVHam3 = (2.3+0.1im)*StepHam 
     @test MulVHam3.ts === StepHam.ts
     @test MulVHam3.ts === Ham.ts
 
     for (i,f) in enumerate(Ham.fs)
         @test MulVHam3.fvals[i] == (2.3+0.1im)*f(t)
     end   

     @test ishermitian(MulVHam3) == false
     @test isskewhermitian(MulVHam3) == false



end


@testset "is_time_dependent" begin
    nsites = 3
    atoms = [(1, i) for i = 1:nsites]
    params = [nothing, 1.0, [0.1 for _ in 1:nsites], sin, [sin for _ in 1:nsites]]
    @testset "is_time_dependent (2-level)" begin
        @testset "Ω=$Ω, Δ=$Δ, ϕ=$ϕ" for Ω in params[2:end], Δ in params, ϕ in params
            h = rydberg_h(atoms; Ω, Δ, ϕ)
            @test !(BloqadeExpr.is_time_dependent(h)) == all(x->(BloqadeExpr.is_const_param(x) || isnothing(x)), [Ω, Δ, ϕ])
        end
    end;
    @testset "is_time_dependent (3-level)" begin
        @testset "Ω_hf=$Ω_hf, Δ_hf=$Δ_hf, ϕ_hf=$ϕ_hf" for Ω_hf in params[2:2], Δ_hf in params[3:4], ϕ_hf in params[5:5]
            @testset "Ω_r=$Ω_r, Δ_r=$Δ_r, ϕ_r=$ϕ_r" for Ω_r in params[2:end], Δ_r in params, ϕ_r in params
                h = rydberg_h_3(atoms; Ω_hf, Δ_hf, ϕ_hf, Ω_r, Δ_r, ϕ_r)
                @test !(BloqadeExpr.is_time_dependent(h)) == all(x->(BloqadeExpr.is_const_param(x) || isnothing(x)), [Ω_hf, Δ_hf, ϕ_hf, Ω_r, Δ_r, ϕ_r])
            end
        end
    end;
    @test !BloqadeExpr.is_time_dependent(nothing)
    @test !BloqadeExpr.is_time_function(SumOfN(nsites))
end
