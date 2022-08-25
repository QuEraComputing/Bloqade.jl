using BloqadeQMC
using LinearAlgebra
using JSON

using CSV
using Tables

function make_hamiltonian_matrix(H::Union{AbstractLTFIM, AbstractTFIM})
    if !haslongitudinalfield(H)
        J, hx, hz = Symmetric(H.J), H.hx, zero(H.hx)
    else
        J, hx, hz = Symmetric(H.J), H.hx, H.hz
    end

    N = nspins(H)
    dim = 2^N

    Hamiltonian = zeros(dim,dim)   #This is your 2D Hamiltonian matrix

    for Ket = 0:dim-1  #Loop over Hilbert Space
        Diagonal = 0.0
        # interaction piece
        for SpinIndex = 0:N-1  #Loop over spin index (zero indexing)
            Spin1 = 2*((Ket>>SpinIndex)&1) - 1

            for NextIndex = SpinIndex+1:N-1
                Spin2 = 2*((Ket>>NextIndex)&1) - 1
                Diagonal += J[SpinIndex+1, NextIndex+1]*Spin1*Spin2
            end

            if hz isa AbstractArray
                Diagonal -= hz[SpinIndex+1]*Spin1
            else
                Diagonal -= hz*Spin1
            end

            bit = 2^SpinIndex   #The "label" of the bit to be flipped
            Bra = Ket ⊻ bit    #Binary XOR flips the bit
            Hamiltonian[Bra+1,Ket+1] = -hx[SpinIndex+1]
        end

        Hamiltonian[Ket+1,Ket+1] = Diagonal
    end

    return Hermitian(Hamiltonian)
end

function make_hamiltonian_matrix(H::AbstractRydberg)
    V, Ω, δ = H.V, H.Ω, H.δ

    N = nspins(H)
    dim = 2^N

    Hamiltonian = zeros(dim,dim)   #This is your 2D Hamiltonian matrix

    for Ket = 0:dim-1  #Loop over Hilbert Space
        Diagonal = 0.0
        # interaction piece
        for SpinIndex = 0:N-2  #Loop over spin index (zero indexing)
            Spin1 = ((Ket>>SpinIndex)&1)

            for NextIndex = SpinIndex+1:N-1
                Spin2 = ((Ket>>NextIndex)&1)
                Diagonal += V[SpinIndex+1, NextIndex+1]*Spin1*Spin2
            end

            if δ isa AbstractArray
                Diagonal -= δ[SpinIndex+1]*Spin1
            else
                Diagonal -= δ*Spin1
            end

            bit = 2^SpinIndex   #The "label" of the bit to be flipped
            Bra = Ket ⊻ bit     #Binary XOR flips the bit
            Hamiltonian[Bra+1,Ket+1] = -0.5*Ω[SpinIndex+1]
        end

        # did not loop over Nth spin in SpinIndex loop
        SpinN = ((Ket>>(N-1))&1)
        if δ isa AbstractArray
            Diagonal -= δ[N]*SpinN
        else
            Diagonal -= δ*SpinN
        end

        bit = 2^(N-1)   #The "label" of the bit to be flipped
        Bra = Ket ⊻ bit     #Binary XOR flips the bit
        Hamiltonian[Bra+1,Ket+1] = -0.5*Ω[N]
        Hamiltonian[Ket+1,Ket+1] = Diagonal
    end

    return Hermitian(Hamiltonian)
end

R_b = 1.7
delta = 3.3
Ω = 1.0
n1 = 2
n2 = 2
t = 1.0


lat = Kagome(t, n1, n2, (false, true); trunc=Inf)
@show lat.distance_matrix
CSV.write("BC=(false,true).csv", Tables.table(lat.distance_matrix), writeheader=false)
println()
lat = Kagome(t, n1, n2, (true, false); trunc=Inf)
CSV.write("BC=(true,false).csv", Tables.table(lat.distance_matrix), writeheader=false)
@show lat.distance_matrix
println()
lat = Kagome(t, n1, n2, (true, true); trunc=Inf)
CSV.write("BC=(true,true).csv", Tables.table(lat.distance_matrix), writeheader=false)
@show lat.distance_matrix
println()

@show R_b, delta
H = Rydberg(lat, R_b, Ω, delta)
N = nspins(H)

Dim = 2^N

SumSz = dropdims(sum(@. (2 * (((0:Dim-1) >> (0:N-1)') & 1) - 1); dims=2); dims=2)
H_mat = make_hamiltonian_matrix(H)
diag = eigen(H_mat)

magnetization = SumSz' * abs2.(diag.vectors)

mag_density = magnetization[1] / N
E_density = diag.values[1] / N
@show E_density
@show mag_density

d = Dict{Symbol, Float64}(:energy_density => E_density, :mag_density => mag_density)
open("../../qmc_data/kagome_nematic_staggered/exact_energy_magnetization_Rb=$(R_b)_delta=$(delta).json", "w") do f
    JSON.print(f, Dict([k => d[k] for k in keys(d)]), 2)
end

println()

#=
for R_b in R_bs
    for delta in δ
        @show R_b, delta
        H = Rydberg(lat, R_b, Ω, delta)
        N = nspins(H)

        Dim = 2^N

        SumSz = dropdims(sum(@. (2 * (((0:Dim-1) >> (0:N-1)') & 1) - 1); dims=2); dims=2)
        H_mat = make_hamiltonian_matrix(H)
        diag = eigen(H_mat)
        
        magnetization = SumSz' * abs2.(diag.vectors)

        mag_density = magnetization[1] / N
        E_density = diag.values[1] / N

        d = Dict{Symbol, Float64}(:energy_density => E_density, :mag_density => mag_density)
        open("../../qmc_data/kagome_nematic_staggered/exact_energy_magnetization_Rb=$(R_b)_delta=$(delta).json", "w") do f
            JSON.print(f, Dict([k => d[k] for k in keys(d)]), 2)
        end

        println()

    end
end
=#
