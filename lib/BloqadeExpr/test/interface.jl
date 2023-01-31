using Test
using BloqadeExpr
using YaoBlocks
using YaoBlocks.Optimise
using BloqadeLattices


@testset "process_atom_positions" begin
    position = [(1, 2), (2, 3)]
    lattice = parallelepiped_region(SquareLattice(),(0,2),(2,0);scale=5)
    @test BloqadeExpr.process_atom_positions(position) == position
    @test BloqadeExpr.process_atom_positions(lattice) == lattice
    
end

@testset "rydberg_h" begin
    atom_list = [(1, 2), (2, 3)]
    lattice = parallelepiped_region(ChainLattice(),(2,);scale=5)

    for positions in [atom_list,lattice]
        h = RydInteract(; atoms = positions) + SumOfX(; nsites = 2, Ω = 0.5)
        @test BloqadeExpr.add_terms(rydberg_h(positions; Ω = 1.0)) == Optimise.simplify(h)

        h = RydInteract(; atoms = positions) + SumOfX(; nsites = 2, Ω = 0.5) - SumOfN(; nsites = 2, Δ = 0.2)
        @test BloqadeExpr.add_terms(rydberg_h(positions; Ω = 1.0, Δ = 0.2)) == Optimise.simplify(h)

        h = RydInteract(; atoms = positions) + SumOfXPhase(; nsites = 2, Ω = 0.5, ϕ = 0.1) - SumOfN(; nsites = 2, Δ = 0.2)
        @test BloqadeExpr.add_terms(rydberg_h(positions; Ω = 1.0, ϕ = 0.1, Δ = 0.2)) == Optimise.simplify(h)

        h = RydInteract(; atoms = positions) + SumOfX(; nsites = 2, Ω = [2.0, 1.0])
        @test BloqadeExpr.add_terms(rydberg_h(positions; Ω = [4.0, 2.0])) == Optimise.simplify(h)
    end

end

@testset "rydberg_h_3" begin
    atom = [(0.0, 0.0)]
    params = [0.0, pi/2, nothing]
    n_to_0(x) = isnothing(x) ? 0.0 : x
    for Ω_hf in params, ϕ_hf in params, Δ_hf in params, Ω_r in params, ϕ_r in params, Δ_r in params
        h = rydberg_h_3(atom; Ω_hf, ϕ_hf, Δ_hf, Ω_r, ϕ_r, Δ_r)
        m = zeros(ComplexF64, 3, 3)
        Ω_hf, ϕ_hf, Δ_hf, Ω_r, ϕ_r, Δ_r = n_to_0.([Ω_hf, ϕ_hf, Δ_hf, Ω_r, ϕ_r, Δ_r])
        m[1, 2] = Ω_hf/2*exp(im*ϕ_hf)
        m[2, 1] = Ω_hf/2*exp(-im*ϕ_hf)
        m[2, 2] = -Δ_hf
        m[2, 3] = Ω_r/2*exp(im*ϕ_r)
        m[3, 2] = Ω_r/2*exp(-im*ϕ_r)
        m[3, 3] = -Δ_hf - Δ_r
        @test BloqadeExpr.mat(h) ≈ m
    end
end

@testset "attime" begin
    positions = [(1, 2), (2, 3)]
    h1 = RydInteract(; atoms = positions) + SumOfX(; nsites = 2, Ω = sin)
    h2 = RydInteract(; atoms = positions) + SumOfX(; nsites = 2, Ω = sin(0.1))
    @test isreal(h1)
    @test isreal(h2)
    @test h1 |> attime(0.1) == h2

    h1 = RydInteract(; atoms = positions) + SumOfX(; nsites = 2, Ω = sin) - SumOfN(; nsites = 2, Δ = cos)
    h2 = RydInteract(; atoms = positions) + SumOfX(; nsites = 2, Ω = sin(0.1)) - SumOfN(; nsites = 2, Δ = cos(0.1))
    @test h1 |> attime(0.1) == h2

    h1 = RydInteract(; atoms = positions) + SumOfXPhase(; nsites = 2, Ω = 0.5, ϕ = cos) - SumOfN(; nsites = 2, Δ = 0.2)
    h2 =
        RydInteract(; atoms = positions) + SumOfXPhase(; nsites = 2, Ω = 0.5, ϕ = cos(0.1)) -
        SumOfN(; nsites = 2, Δ = 0.2)
    @test !isreal(h1)
    @test h1 |> attime(0.1) == h2

    h1 = RydInteract(; atoms = positions) + SumOfX(; nsites = 2, Ω = [sin, cos])
    h2 = RydInteract(; atoms = positions) + SumOfX(; nsites = 2, Ω = [sin(0.1), cos(0.1)])
    @test h1 |> attime(0.1) == h2
end


@testset "get_rydberg_params" begin
    atoms = [(i,i) for i in 1:10]

    values = [nothing,1,t->t^2,rand(10),[t->rand()*t for i in 1:10]]

    for ϕ in values, Ω in values, Δ in values
        h = rydberg_h(atoms,ϕ=ϕ,Ω=Ω,Δ=Δ)
        # catches the this weird edge case 
        Ω = (isnothing(Ω) && !isnothing(ϕ) ? 0 : Ω)
        @test (atoms=atoms,ϕ=ϕ,Ω=Ω,Δ=Δ) == get_rydberg_params(h)
    end
end

