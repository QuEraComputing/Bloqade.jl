using Test
using Adapt
using Random
using RydbergEmulator
using Graphs: SimpleGraph, add_edge!
using SparseArrays
using OrderedCollections
using LuxurySparse
using BitBasis
using Printf
using Unitful
using Yao
using Unitful: μm, mm, μs, ns, MHz, GHz
using Yao.ConstGate: P0, P1
using RydbergEmulator: Negative, PrecisionAdaptor, sparse_skeleton_csc

Random.seed!(1234)

function rydinteract(atoms, C)
    n = length(atoms)
    terms = []
    for i in 1:n, j in 1:n
        if i != j
            push!(terms, C/(2 * RydbergEmulator.distance(atoms[i], atoms[j])^6) * kron(n, i=>P1, j=>P1))
        end
    end
    return sum(terms)
end

if !isdefined(@__MODULE__, :test_graph)
    include("utils.jl")
end

@testset "simple graph hamiltonian subspace" begin
    subspace = blockade_subspace(test_graph)
    show(stdout, MIME"text/plain"(), subspace)
    show(IOContext(stdout, :limit=>true), MIME"text/plain"(), subspace)

    @test collect(keys(subspace)) == sort(test_subspace_v)
    @test collect(values(subspace)) == [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]

    # test hamiltonian creation
    Ω = Float64[2.5, 3.4, 0.2, 1.7, 4.3]
    Δ = Float64[1.2, 3.4, 2.6, 0.2, 1.8]
    ϕ = Float64[0.5, 0.4, -0.2, -1.2, 10.2]

    target = create_test_hamiltonian(Δ, Ω, ϕ)
    h = XTerm(Ω, ϕ) + ZTerm(Δ)
    @test SparseMatrixCSC(h, subspace) ≈ target

    H = SparseMatrixCSC(h, subspace)
    @test H ≈ update_term!(copy(H), h, subspace)
end

@testset "sparse_skeleton_csc($term)" for (op, term) in [Yao.X=>XTerm, P1=>NTerm, Yao.Z=>ZTerm]
    @testset "natoms=$natoms" for natoms in 3:5
        H = mat(sum([1.0 * kron(natoms, k=>op) for k in 1:natoms]))
        if term === NTerm
            # NOTE: we force NTerm matrix to be fully diagonal
            # since only the first entry is zero
            H[1, 1] = 1.0
        end
        H = SparseMatrixCSC(H)
        h = term(natoms, 1.0)

        @testset "fullspace" begin
            colptr, rowval = sparse_skeleton_csc(h)
            @test H.colptr == colptr
            @test H.rowval == rowval
        end

        @testset "subspace trial=$trial" for trial in 1:3
            space = Subspace(randperm(1<<natoms)[1:1<<(natoms-2)].-1)
            colptr, rowval = sparse_skeleton_csc(h, space)

            M = H[vec(space).+1, vec(space).+1]
            @test M.colptr == colptr
            @test M.rowval == rowval
        end
    end
end

@testset "SparseMatrixCSC($term)" for ((a, term), (b, op)) in [
        (4.0=>XTerm, 2.0=>Yao.X),
        (2.0=>ZTerm, 2.0=>Yao.Z),
        (2.0=>NTerm, 2.0=>P1),
    ]

    @testset "natoms=$natoms" for natoms in 3:5
        H = mat(sum([b * kron(natoms, k=>op) for k in 1:natoms]))
        h = term(natoms, a)
        @testset "fullspace" begin
            @test SparseMatrixCSC(h) ≈ H
            @testset "Tv=$Tv" for Tv in [Float32, Float64, ComplexF32, ComplexF64]
                @test SparseMatrixCSC{Tv}(h) ≈ H

                @testset "Ti=$Ti" for Ti in [Cint, Int]
                    @test SparseMatrixCSC{Tv, Ti}(h) ≈ H
                end
            end
        end # fullspace

        @testset "subspace trial=$trial" for trial in 1:3
            space = Subspace(randperm(1<<natoms)[1:1<<(natoms-2)].-1)
            M = H[vec(space).+1, vec(space).+1]
            @test SparseMatrixCSC(h, space) ≈ M
            @testset "Tv=$Tv" for Tv in [Float32, Float64, ComplexF32, ComplexF64]
                @test SparseMatrixCSC{Tv}(h, space) ≈ M

                @testset "Ti=$Ti" for Ti in [Cint, Int]
                    @test SparseMatrixCSC{Tv, Ti}(h, space) ≈ M
                end
            end
        end # subspace trial
    end
end

@testset "SparseMatrixCSC(RydInteract) natoms=$natoms" for natoms in 3:5
    atoms = square_lattice(natoms, 0.8)
    H = mat(rydinteract(atoms, 1.0))
    h = RydInteract(atoms, 1.0)
    @testset "fullspace" begin
        @test SparseMatrixCSC(h) ≈ H
        @testset "Tv=$Tv" for Tv in [Float32, Float64, ComplexF32, ComplexF64]
            @test SparseMatrixCSC{Tv}(h) ≈ H
            @testset "Ti=$Ti" for Ti in [Cint, Int]
                @test SparseMatrixCSC{Tv, Ti}(h) ≈ H
            end
        end
    end # fullspace

    @testset "subspace trial=$trial" for trial in 1:3
        space = Subspace(randperm(1<<natoms)[1:1<<(natoms-2)].-1)
        M = H[vec(space).+1, vec(space).+1]
        @test SparseMatrixCSC(h, space) ≈ M
        @testset "Tv=$Tv" for Tv in [Float32, Float64, ComplexF32, ComplexF64]
            @test SparseMatrixCSC{Tv}(h, space) ≈ M

            @testset "Ti=$Ti" for Ti in [Cint, Int]
                @test SparseMatrixCSC{Tv, Ti}(h, space) ≈ M
            end
        end
    end # subspace trial
end

@testset "X term" begin
    H = mat(sum([2.0 * kron(5, k=>Yao.X) for k in 1:5]))
    h = XTerm(5, 4.0)
    @test h(0.1) == h
    @test SparseMatrixCSC(h) ≈ H
    @test update_term!(SparseMatrixCSC(h), h) ≈ H
    @test update_term!(SparseMatrixCSC(h), h, FullSpace()) ≈ H
    display(h)

    Ωs = rand(5)
    h = XTerm(Ωs)
    H = mat(sum([Ω/2 * kron(5, k=>Yao.X) for (k, Ω) in enumerate(Ωs)]))
    @test SparseMatrixCSC(h) ≈ H
    @test update_term!(SparseMatrixCSC(h), h) ≈ H
    @test update_term!(SparseMatrixCSC(h), h, FullSpace()) ≈ H
    @test eltype(h) == Float64
    display(h)

    Ωs = rand(5)
    ϕs = rand(5)
    h = XTerm(Ωs, ϕs)
    display(h)

    ϕs = rand(5)
    h = XTerm(2.02, ϕs)
    display(h)

    h = XTerm(5, sin)
    @test isreal(h) == true
    @test h(0.1) isa XTerm
    @test h(0.1).Ωs == sin(0.1)
    @test h(0.1).ϕs === nothing

    h = XTerm(5, sin, cos)
    @test h(0.1) isa XTerm
    @test h(0.1).Ωs == sin(0.1)
    @test h(0.1).ϕs == cos(0.1)

    @test isreal(h) == false

    display(h)

    @test_throws AssertionError XTerm(5, [1, 2, 3, 4])
end

@testset "Z term" begin
    H = mat(sum([2.0 * kron(5, k=>Yao.Z) for k in 1:5]))
    h = ZTerm(5, 2.0)
    @test SparseMatrixCSC(h) ≈ H
    @test update_term!(SparseMatrixCSC(h), h) ≈ H

    display(h)

    Δs = rand(5)
    h = ZTerm(Δs)
    H = mat(sum([Δs[k] * kron(5, k=>Yao.Z) for k in 1:5]))
    @test SparseMatrixCSC(h) ≈ H
    @test update_term!(SparseMatrixCSC(h), h) ≈ H
    @test eltype(h) == Float64
    display(h)

    h = ZTerm(5, sin)
    @test isreal(h) == true
    @test h(0.1).Δs == sin(0.1)

    h = ZTerm(3, [sin, cos, tanh])
    @test h(0.1).Δs[1] == sin(0.1)
    @test h(0.1).Δs[2] == cos(0.1)
    @test h(0.1).Δs[3] == tanh(0.1)
end

@testset "N term" begin
    H = mat(sum([2.0 * kron(5, k=>P1) for k in 1:5]))
    h = NTerm(5, 2.0)
    @test SparseMatrixCSC(h) ≈ H
    @test update_term!(SparseMatrixCSC(h), h) ≈ H

    display(h)

    Δs = rand(5)
    h = NTerm(Δs)
    H = mat(sum([Δs[k] * kron(5, k=>P1) for k in 1:5]))
    @test SparseMatrixCSC(h) ≈ H
    @test update_term!(SparseMatrixCSC(h), h) ≈ H
    @test eltype(h) == Float64
    display(h)

    h = NTerm(5, sin)
    @test isreal(h) == true
    @test h(0.1).Δs == sin(0.1)

    h = NTerm(3, [sin, cos, tanh])
    @test h(0.1).Δs[1] == sin(0.1)
    @test h(0.1).Δs[2] == cos(0.1)
    @test h(0.1).Δs[3] == tanh(0.1)
end


@testset "rydberg interact term" begin
    atoms = RydbergEmulator.square_lattice(4, 0.8)
    H = mat(rydinteract(atoms, 2.0))
    h = RydInteract(atoms, 2.0)
    @test SparseMatrixCSC(h) ≈ H
    @test update_term!(SparseMatrixCSC(h), h) ≈ H
    @test eltype(h) == Float64

    display(h)
end

@testset "composite term" begin
    atoms = RydbergEmulator.square_lattice(5, 0.8)
    h = XTerm(5, 4.0) + RydInteract(atoms, 2.0) + ZTerm(5, 1.0)
    H = rydinteract(atoms, 2.0) +
        sum([1.0 * kron(5, k=>Yao.Z) for k in 1:5]) +
        sum([4.0/2 * kron(5, k=>Yao.X) for k in 1:5])
    H = mat(H)
    @test SparseMatrixCSC(h) ≈ H
    @test update_term!(SparseMatrixCSC(h), h) ≈ H

    display(h)

    h1 = XTerm(5, 4.0) + RydInteract(atoms, 2.0)
    @test isreal(h1) == true
    @test ZTerm(5, 1.0) + h1 == Hamiltonian((ZTerm(5, 1.0), h1.terms...))
    @test h1 + ZTerm(5, 1.0) == Hamiltonian((h1.terms..., ZTerm(5, 1.0)))
    @test h1 + h1 == Hamiltonian((h1.terms..., h1.terms...))

    h2 = RydInteract(atoms, 2.0) - XTerm(5, 4.0)
    h3 = RydInteract(atoms, 2.0) + XTerm(5, -4.0)
    @test SparseMatrixCSC(h2) == SparseMatrixCSC(h3)

    h4 = RydInteract(atoms, 2.0) + XTerm(5, 4.0) - ZTerm(5,  1.2)
    h5 = RydInteract(atoms, 2.0) + XTerm(5, 4.0) + ZTerm(5, -1.2)
    @test isreal(h4) == true
    @test SparseMatrixCSC(h4) == SparseMatrixCSC(h5)

    h6 = RydInteract(atoms, 2.0) - XTerm([ 1,  2,  3,  4,  5]) - ZTerm(5,  1.2)
    h7 = RydInteract(atoms, 2.0) + XTerm([-1, -2, -3, -4, -5]) + ZTerm(5, -1.2)
    @test SparseMatrixCSC(h6) == SparseMatrixCSC(h7)

    h8 = RydInteract(atoms, 2.0) - XTerm(5, sin) - ZTerm(5,  cos)
    h9 = RydInteract(atoms, 2.0) + Negative(XTerm(5, sin)) + Negative(ZTerm(5,  cos))
    @test h8(1.0) == h9(1.0)

    h = XTerm(5, sin) + NTerm(5, cos) - (XTerm(5, sin) + NTerm(5, cos))
    @test iszero(SparseMatrixCSC(h(1.0)))

    h = NTerm(5, cos) - (XTerm(5, sin) + NTerm(5, cos))
    @test SparseMatrixCSC(h(1.0)) ≈ SparseMatrixCSC(-XTerm(5, sin(1.0)))
end

@testset "XTerm subspace" begin
    Ω = Float64[2.5, 3.4, 0.2, 1.7, 4.3]
    ϕ = Float64[0.5, 0.4, -0.2, -1.2, 10.2]
    h = XTerm(Ω, ϕ)

    subspace = blockade_subspace(test_graph)
    M = SparseMatrixCSC(h, subspace)
    @test update_term!(copy(M), h, subspace) ≈ M
end

@testset "$term subspace" for (term, op) in [(ZTerm, Z), (NTerm, P1)]
    atoms = square_lattice(5, 0.8)
    space = blockade_subspace(atoms, 1.0)
    N = length(space)
    h = term(5, 1.0)
    H = Matrix(sum(kron(5, k=>op) for k in 1:5))
    target_H = H[space.subspace_v .+ 1, space.subspace_v .+ 1]
    test_H = SparseMatrixCSC(h, space)
    @test test_H ≈ target_H
end

@testset "redberg interact term subspace" begin
    atoms = RydbergEmulator.square_lattice(4, 0.8)
    s = blockade_subspace(atoms, 1.5)
    H = mat(rydinteract(atoms, 2.0));
    h = RydInteract(atoms, 2.0)

    @test SparseMatrixCSC(h)[s.subspace_v.+1, s.subspace_v.+1] ≈ SparseMatrixCSC(h, s)
    @test update_term!(SparseMatrixCSC(h, s), h, s) ≈ SparseMatrixCSC(h, s)
end

@testset "utils" begin
    @test RydbergEmulator.to_tuple((1, 2, 3)) == (1, 2, 3)
    @test RydbergEmulator.to_tuple([1, 2, 3]) == (1, 2, 3)
    @test RydbergEmulator.getscalarmaybe([1, 2, 3], 2) == 2
    @test RydbergEmulator.getscalarmaybe(3.21, 5) == 3.21
    @test RydbergEmulator.getscalarmaybe(nothing, 5) == 0
    atoms = rand_atoms(5, 0.3)
    h1 = rydberg_h(atoms, 1.0, 2.21, 3.32, 2.22)
    h2 = RydInteract(atoms, 1.0) + XTerm(5, 2.21, 3.32) - NTerm(5, 2.22)
    @test h1 == h2
end

@testset "units" begin
    atom = RydAtom(1mm, 2μm)
    @test atom[1] == 1000
    @test atom[2] == 2

    h = RydInteract([RydAtom(1.0mm, 2.0μm), RydAtom(1.1mm, 2.2μm)], 2GHz * μm^6)
    @test h.C == 2000
    @test h.atoms[1][1] == 1000
    @test h.atoms[2][1] == 1100

    h = RydInteract([RydAtom(1.0mm, 2.0μm), RydAtom(1.1mm, 2.2μm)])
    @test h.C == 2π * 109.133

    h = XTerm(5, 1.0GHz)
    @test h.Ωs == 1000
    h = XTerm([1.0GHz, 2.0MHz])
    @test h.Ωs[1] == 1000
    @test h.Ωs[2] == 2

    h = XTerm(5, 1.0GHz, 2.0MHz * μs)
    @test h.ϕs == 2.0

    h = ZTerm(5, 1.0GHz)
    @test h.Δs == 1000

    h = ZTerm([1.0GHz, 2.0MHz])
    @test h.Δs[1] == 1000
    @test h.Δs[2] == 2.0
end

check_print(term, line::String) = check_print(term, [line])

function check_print(term, lines::Vector{String})
    buf = IOBuffer()
    io = IOContext(buf, :limit=>true)
    show(io, MIME"text/plain"(), term)
    s = String(take!(buf))

    pass = true
    for (i, l) in enumerate(lines)
        if i != lastindex(lines)
            test_string = l * '\n'
        else
            test_string = l
        end
        pass = pass && occursin(test_string, s)
    end
    return pass
end

@testset "limit printing" begin
    @test check_print(XTerm(10, sin), ["XTerm", " ∑(n=1:10) sin(t)/2 σ^x"])
    @test check_print(XTerm([1, 2, 3, 4, 5, 6, 7, 8, 9, 10]), [
        "XTerm",
        " 1/2 σ^x +",
        " σ^x +",
        " 3/2 σ^x +",
        "   ⋯",
        " 8/2 σ^x +",
        " 9/2 σ^x +",
        " 10/2 σ^x"
    ])

    @test check_print(NTerm(10, sin), ["NTerm", " ∑(n=1:10) sin(t) n"])
end

@testset "adapt" begin
    atoms = square_lattice(10, 0.8)
    h = rydberg_h(atoms, 1.2, 1.3, 2.1)

    new_h = adapt(PrecisionAdaptor(Float32), h)
    @test new_h[1].C isa Float32
    @test new_h[2].Ωs isa Float32
    @test new_h[2].ϕs isa Float32
    @test new_h[3].term.Δs isa Float32
end
