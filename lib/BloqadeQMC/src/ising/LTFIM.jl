using Base.Iterators


abstract type AbstractLTFIM{O <: AbstractOperatorSampler} <: AbstractIsing{O} end

struct GeneralLTFIM{O,M <: UpperTriangular{Float64},Vx <: AbstractVector{Float64},Vz <: AbstractVector{Float64}} <: AbstractLTFIM{O}
    op_sampler::O
    J::M
    hx::Vx
    hz::Vz
    Ns::Int
    Nb::Int
    energy_shift::Float64
end


struct LTFIM{O} <: AbstractLTFIM{O}
    op_sampler::O
    J::Float64
    hx::Float64
    hz::Float64
    Ns::Int
    Nb::Int
    energy_shift::Float64
end
###############################################################################

# LTFIM ops:
#  (-2,i,i) is an off-diagonal site operator h(sigma^+_i + sigma^-_i)
#  (-1,i,i) is a diagonal site operator h
#  (0,0,0) is the identity operator I - NOT USED IN THE PROJECTOR CASE
#  (t,i,j) is a diagonal bond operator J(sigma^z_i sigma^z_j) + hzb(sigma^z_i + sigma^z_j)
#    t denotes the spin config at sites i,j, subtract 1 then convert to binary
#    t = 1 -> 00 -> down-down
#    t = 2 -> 01 -> down-up
#    t = 3 -> 10 -> up-down
#    t = 4 -> 11 -> up-up
#    spin_config(t) = divrem(t - 1, 2)
@inline getbondtype(::AbstractLTFIM, s1::Bool, s2::Bool) = (s1<<1 | s2) + 1
@inline spin_config(::AbstractLTFIM, t::Int)::NTuple{2,Int} = divrem(t - 1, 2)
@inline spin_config(H::AbstractLTFIM, op::NTuple{ISING_OP_SIZE, Int}) = @inbounds spin_config(H, getoperatortype(H, op))

###############################################################################

function make_prob_vector(J::UpperTriangular{T}, hx::AbstractVector{T}, hz::AbstractVector{T}; epsilon=0.0) where T
    @assert length(hx) == length(hz) == size(J, 1) == size(J, 2)

    ops = Vector{NTuple{ISING_OP_SIZE, Int}}()
    p = Vector{T}()
    energy_shift = zero(T)

    for i in eachindex(hx)
        if !iszero(hx[i])
            push!(ops, makediagonalsiteop(AbstractLTFIM, i))
            push!(p, hx[i])
            energy_shift += hx[i]
        end
    end

    Ns = length(hx)
    bond_spins = Set{NTuple{2,Int}}()
    coordination_numbers = zeros(Int, Ns)
    for j in axes(J, 2), i in axes(J, 1)
        if i < j && !iszero(J[i, j])
            push!(bond_spins, (i, j))
            coordination_numbers[i] += 1
            coordination_numbers[j] += 1
        end
    end

    Z = diagonaloperator(AbstractLTFIM)  # since 0 maps to spin down
    I = Diagonal(LinearAlgebra.I, 2)

    # add fictitious bonds if there's a z-field on an "unbonded" site
    while any(i -> iszero(coordination_numbers[i]) && !iszero(hz[i]), 1:Ns)
        # TODO: test this
        f = findfirst(i -> iszero(coordination_numbers[i]) && !iszero(hz[i]), 1:Ns)
        l = findlast(i -> iszero(coordination_numbers[i]) && !iszero(hz[i]), 1:Ns)

        if f == l
            # if there's only one unbonded site,
            # pick some other site to bond it to
            i, _ = first(bond_spins)
            site1, site2 = (f < i) ? (f, i) : (i, f)
        else
            site1, site2 = (f < l) ? (f, l) : (l, f)
        end

        push!(bond_spins, (site1, site2))
        coordination_numbers[site1] += 1
        coordination_numbers[site2] += 1
    end

    for (site1, site2) in bond_spins
        # by this point we can assume site1 < site2
        hzb1 = hz[site1] / coordination_numbers[site1]
        hzb2 = hz[site2] / coordination_numbers[site2]
        local_H = J[site1, site2]*kron(Z, Z) - hzb1*kron(Z, I) - hzb2*kron(I, Z)

        p_spins = -diag(local_H)
        C = abs(min(0, minimum(p_spins))) + epsilon
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

function LTFIM(dims::NTuple{N, Int}, J::Float64, hx::Float64, hz::Float64, pbc=true) where N
    bond_spins, Ns, Nb = lattice_bond_spins(dims, pbc)
    J_, hx_ = make_uniform_tfim(bond_spins, Ns, J, hx)
    ops, p, energy_shift = make_prob_vector(J_, hx_, hz*ones(Ns))
    op_sampler = ImprovedOperatorSampler(AbstractLTFIM, ops, p)
    return LTFIM{typeof(op_sampler)}(op_sampler, J, hx, hz, Ns, Nb, energy_shift)
end

function GeneralLTFIM(J::UpperTriangular{Float64}, hx::AbstractVector{Float64}, hz::AbstractVector{Float64}) where N
    ops, p, energy_shift = make_prob_vector(J, hx, hz)
    op_sampler = ImprovedOperatorSampler(AbstractLTFIM, ops, p)
    return GeneralLTFIM{typeof(op_sampler),typeof(J),typeof(hx),typeof(hz)}(op_sampler, J, hx, hz, Ns, Nb, energy_shift)
end

function GeneralLTFIM(dims::NTuple{N, Int}, J::Float64, hx::Float64, hz::Float64, pbc=true) where N
    bond_spins, Ns, Nb = lattice_bond_spins(dims, pbc)
    J_, hx_ = make_uniform_tfim(bond_spins, Ns, J, hx)
    return GeneralLTFIM(J_, hx_, hz*ones(Ns))
end

total_hx(H::GeneralLTFIM)::Float64 = sum(H.hx)
total_hx(H::LTFIM)::Float64 = H.hx * nspins(H)
haslongitudinalfield(H::AbstractLTFIM) = !iszero(H.hz)