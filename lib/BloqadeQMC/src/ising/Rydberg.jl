using Base.Iterators
using BloqadeLattices: rydberg_interaction_matrix, BoundedLattice

abstract type AbstractRydberg{O <: AbstractOperatorSampler} <: AbstractLTFIM{O} end

struct Rydberg{O,M <: UpperTriangular,UΩ <: AbstractVector{Float64}, Uδ <: AbstractVector{Float64}, A} <: AbstractRydberg{O}
    op_sampler::O
    V::M          # interaction matrix
    Ω::UΩ 
    δ::Uδ
    atoms::A          # hoping to include typing ::Union{Vector,BoundedLattice} during refactoring. This actually does break stuff!
    energy_shift::Float64
end

nspins(H::Rydberg) = length(H.atoms)

@inline diagonaloperator(::Type{<:AbstractRydberg}) = Diagonal([0, 1])
@inline diagonaloperator(H::AbstractRydberg) = diagonaloperator(typeof(H))


function make_prob_vector(H::Type{<:AbstractRydberg}, V::UpperTriangular, Ω::AbstractVector{T}, δ::AbstractVector{T}; epsilon=0.0) where T
    @assert length(Ω) == length(δ) == size(V, 1) == size(V, 2)
    @assert (0.0 <= epsilon <= 1.0) "epsilon must be in the range [0, 1]!"

    ops = Vector{NTuple{ISING_OP_SIZE, Int}}()
    p = Vector{T}()
    energy_shift = zero(T)

    for i in eachindex(Ω)
        if !iszero(Ω[i])
            push!(ops, makediagonalsiteop(AbstractLTFIM, i))
            push!(p, Ω[i] / 2)
            energy_shift += Ω[i] / 2
        end
    end

    Ns = length(Ω)
    bond_spins = Set{NTuple{2,Int}}()
    coordination_numbers = zeros(Int, Ns)
    for j in axes(V, 2), i in axes(V, 1)
        if i < j && !iszero(V[i, j])
            push!(bond_spins, (i, j))
            coordination_numbers[i] += 1
            coordination_numbers[j] += 1
        end
    end

    n = diagonaloperator(H)
    I = Diagonal(LinearAlgebra.I, 2)

    # TODO: add fictitious bonds if there's a z-field on an "unbonded" site
    for (site1, site2) in bond_spins
        # by this point we can assume site1 < site2
        δb1 = δ[site1] / coordination_numbers[site1]
        δb2 = δ[site2] / coordination_numbers[site2]
        local_H = V[site1, site2]*kron(n, n) - δb1*kron(n, I) - δb2*kron(I, n)

        p_spins = -diag(local_H)
        C = abs(min(0, minimum(p_spins))) + epsilon*abs(minimum(p_spins[2:end]))
        #dont use the zero matrix element for the epsilon shift
        p_spins .+= C
        energy_shift += C

        for (t, p_t) in enumerate(p_spins)
            push!(p, p_t)
            push!(ops, (2, t, length(p), site1, site2))
        end
    end

    return ops, p, energy_shift
end

total_hx(H::Rydberg)::Float64 = sum(H.Ω) / 2
haslongitudinalfield(H::AbstractRydberg) = !iszero(H.δ)

function make_vector(param::AbstractVector, Ns::Int)
    return param
end

function make_vector(param::Real, Ns::Int)
    return param*ones(Ns)
end


function rydberg_qmc(h::RydbergHamiltonian) 
    atoms,ϕ,Ω,Δ = get_rydberg_params(h)

    if !isnothing(ϕ)
        error("SSE QMC currently does not support a non-zero laser phase ϕ")
    end

    if is_time_function(Ω) || is_time_function(Δ)
        error("SSE QMC currently does not support time-dependent waveforms")
    end

    Ns = length(atoms)
    C = h.rydberg_term.C
    V = rydberg_interaction_matrix(atoms, C)

    Ω_N = make_vector(Ω,Ns)
    Δ_N = make_vector(Δ,Ns)

    ops, p, energy_shift = make_prob_vector(AbstractRydberg, V, Ω_N, Δ_N, epsilon=0.0)
    op_sampler = ImprovedOperatorSampler(AbstractLTFIM, ops, p)
    return Rydberg{typeof(op_sampler), typeof(V), typeof(Ω_N), typeof(Δ_N), typeof(atoms)}(op_sampler, V, Ω_N, Δ_N, atoms, energy_shift)
end
