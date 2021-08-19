using Test
using RydbergEmulator
using LightGraphs: SimpleGraph, add_edge!
using SparseArrays
using OrderedCollections
using LuxurySparse
using BitBasis
using Printf
using Unitful
using Yao
using Unitful: μm, mm, μs, ns, MHz, GHz
using Yao.ConstGate: P0, P1

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

    # @testset "cuda" begin
    #     if CUDA.functional()
    #         using CUDA
    #         using CUDA.CUSPARSE
    #         dΩ = cu(Ω)
    #         dϕ = cu(ϕ)
    #         dΔ = cu(Δ)
    #         h = XTerm(dΩ, dϕ) + ZTerm(dΔ)
    #         H = SparseMatrixCSC(h, subspace)
    #         dH = CuSparseMatrixCSR(H)
    #         ds = cu(subspace)
    #         update_term!(dH, cu(h), ds)
    #         @test isapprox(SparseMatrixCSC(dH), H; rtol=1e-7)
    #     end
    # end
end

@testset "X term" begin
    H = mat(sum([2.0 * kron(5, k=>Yao.X) for k in 1:5]))
    h = XTerm(5, 4.0)
    @test SparseMatrixCSC(h) ≈ H
    @test update_term!(SparseMatrixCSC(h), h) ≈ H
    display(h)

    Ωs = rand(5)
    h = XTerm(Ωs)
    H = mat(sum([Ω/2 * kron(5, k=>Yao.X) for (k, Ω) in enumerate(Ωs)]))
    @test SparseMatrixCSC(h) ≈ H
    @test update_term!(SparseMatrixCSC(h), h) ≈ H
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
    @test h(0.1) isa XTerm
    @test h(0.1).Ωs == sin(0.1)
    @test h(0.1).ϕs === nothing

    h = XTerm(5, sin, cos)
    @test h(0.1) isa XTerm
    @test h(0.1).Ωs == sin(0.1)
    @test h(0.1).ϕs == cos(0.1)

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
    @test ZTerm(5, 1.0) + h1 == Hamiltonian((ZTerm(5, 1.0), h1.terms...))
    @test h1 + ZTerm(5, 1.0) == Hamiltonian((h1.terms..., ZTerm(5, 1.0)))
    @test h1 + h1 == Hamiltonian((h1.terms..., h1.terms...))

    h2 = RydInteract(atoms, 2.0) - XTerm(5, 4.0)
    h3 = RydInteract(atoms, 2.0) + XTerm(5, -4.0)
    @test SparseMatrixCSC(h2) == SparseMatrixCSC(h3)

    h4 = RydInteract(atoms, 2.0) + XTerm(5, 4.0) - ZTerm(5,  1.2)
    h5 = RydInteract(atoms, 2.0) + XTerm(5, 4.0) + ZTerm(5, -1.2)
    @test SparseMatrixCSC(h4) == SparseMatrixCSC(h5)

    h6 = RydInteract(atoms, 2.0) - XTerm([ 1,  2,  3,  4,  5]) - ZTerm(5,  1.2)
    h7 = RydInteract(atoms, 2.0) + XTerm([-1, -2, -3, -4, -5]) + ZTerm(5, -1.2)
    @test SparseMatrixCSC(h6) == SparseMatrixCSC(h7)

    h8 = RydInteract(atoms, 2.0) - XTerm(5, sin) - ZTerm(5,  cos)
    h9 = RydInteract(atoms, 2.0) + XTerm(5, x->-sin(x)) + ZTerm(5, x->-cos(x))
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
    H = SparseMatrixCSC(h, subspace)
    @test update_term!(copy(H), h, subspace) ≈ H
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
    h2 = RydInteract(atoms, 1.0) + XTerm(5, 2.21, 3.32) + ZTerm(5, 2.22)
    @test h1 == h2
end

@testset "Yao.mat" begin
    t = simple_rydberg(5, 0.8)
    @test mat(t) ≈ SparseMatrixCSC(t)
    @test mat(t, test_subspace) ≈ SparseMatrixCSC(t, test_subspace)
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
