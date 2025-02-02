using Base.Iterators

export Rydberg

# changed type tree s.t. AbstractRydberg is a direct subtype of Hamiltonian... ACTUALLY NO NEED FOR ABSTRACT AT ALL?
# abstract type AbstractRydberg{O <: AbstractOperatorSampler} <: Hamiltonian{2,O} end

struct Rydberg{O,M <: UpperTriangular{Float64},UΩ <: AbstractVector{Float64}, Uδ <: AbstractVector{Float64}, L <: Lattice} # <: AbstractRydberg{O}
    op_sampler::O
    V::M
    Ω::UΩ
    δ::Uδ
    lattice::L
    energy_shift::Float64
end

nspins(H::Rydberg) = nspins(H.lattice)


function make_prob_vector(H::Type{<:Rydberg}, V::UpperTriangular{T}, Ω::AbstractVector{T}, δ::AbstractVector{T}; epsilon=0.0) where T
    @assert length(Ω) == length(δ) == size(V, 1) == size(V, 2)
    @assert (0.0 <= epsilon <= 1.0) "epsilon must be in the range [0, 1]!"

    ops = Vector{NTuple{OP_SIZE, Int}}()
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


###############################################################################

# function Rydberg(dims::NTuple{D, Int}, R_b, Ω, δ; pbc=true, trunc::Int=0, epsilon::Float64=0) where D
#     if D == 1
#         lat = Chain(dims[1], 1.0, pbc)
#     elseif D == 2
#         lat = Rectangle(dims[1], dims[2], 1.0, 1.0, pbc)
#     else
#         error("Unsupported number of dimensions. 1- and 2-dimensional lattices are supported.")
#     end
#     return Rydberg(lat, R_b, Ω, δ; trunc=trunc, epsilon=epsilon)
# end


# We  only want one Rydberg() interface, integrated into BloqadeExpr. 
# Probably, we'll want one function generate_interaction_matrix() that creates V and then feeds that matrix into Rydberg(). 
# Check with Phil how that function is coming along.

function Rydberg(lat::Lattice, R_b::Float64, Ω::Float64, δ::Float64; trunc::Int=0, epsilon=0)
    println("Hello!")
    Ns = nspins(lat)
    V = zeros(Float64, Ns, Ns)

    if trunc > 0
        _dist = sort!(collect(Set(lat.distance_matrix)))
        uniq_dist = Vector{Float64}(undef, 0)
        for i in eachindex(_dist)
            if length(uniq_dist) == 0
                push!(uniq_dist, _dist[i])
            elseif !(last(uniq_dist) ≈ _dist[i])
                push!(uniq_dist, _dist[i])
            end
        end
        smallest_k = sort!(uniq_dist)[1:(trunc+1)]
        dist = copy(lat.distance_matrix)
        for i in eachindex(dist)
            if dist[i] > last(smallest_k) && !(dist[i] ≈ last(smallest_k))
                dist[i] = zero(dist[i])
            end
        end
    elseif lat isa Rectangle && all(lat.PBC)
        V = zeros(Ns, Ns)
        K = 3  # figure out an efficient way to set this dynamically

        dist = zeros(Ns, Ns)
        for v2 in -K:K, v1 in -K:K
            dV = zeros(Ns, Ns)
            for x2 in axes(dV, 2), x1 in axes(dV, 1)
                i1, j1 = divrem(x1 - 1, lat.n1)
                i2, j2 = divrem(x2 - 1, lat.n1)
                r = [i2 - i1 + v1*lat.n1, j2 - j1 + v2*lat.n2]
                dV[x1, x2] += Ω * (R_b/norm(r, 2))^6
            end
            # @show v2, v1, maximum(abs, dV)

            V += dV
        end

        V = (V + V') / 2  # should already be symmetric but just in case

        return Rydberg(UpperTriangular(triu!(V, 1)), Ω*ones(Ns), δ*ones(Ns), lat; epsilon=epsilon)
    else
        dist = lat.distance_matrix
    end

    @inbounds for i in 1:(Ns-1)
        for j in (i+1):Ns
            # a zero entry in distance_matrix means there should be no bond
            V[i, j] = dist[i, j] != 0.0 ? Ω * (R_b / dist[i, j])^6 : 0.0
        end
    end
    V = UpperTriangular(triu!(V, 1))

    return Rydberg(V, Ω*ones(Ns), δ*ones(Ns), lat; epsilon=epsilon)
end

function Rydberg(V::AbstractMatrix{T}, Ω::AbstractVector{T}, δ::AbstractVector{T}, lattice::Lattice; epsilon=zero(T)) where T
    ops, p, energy_shift = make_prob_vector(Rydberg, V, Ω, δ, epsilon=epsilon)
    op_sampler = ImprovedOperatorSampler(AbstractLTFIM, ops, p)
    return Rydberg{typeof(op_sampler), typeof(V), typeof(Ω), typeof(δ), typeof(lattice)}(op_sampler, V, Ω, δ, lattice, energy_shift)
end

# Check whether the two functions below are needed in updates.

total_hx(H::Rydberg)::Float64 = sum(H.Ω) / 2
haslongitudinalfield(H::Rydberg) = !iszero(H.δ)
