abstract type AbstractIsing{O <: AbstractOperatorSampler} <: Hamiltonian{2,O} end
abstract type AbstractTFIM{O <: AbstractOperatorSampler} <: AbstractIsing{O} end

# H = -J ∑_{⟨ij⟩} σ_i σ_j - h ∑_i σ^x_i
struct NearestNeighbourTFIM{O} <: AbstractTFIM{O}
    op_sampler::O
    J::Float64
    hx::Float64
    Ns::Int
    energy_shift::Float64
end


# H = ∑_{ij} J_ij σ_i σ_j - ∑_i h_i σ^x_i
struct TFIM{O,M <: UpperTriangular{Float64},V <: AbstractVector{Float64}} <: AbstractTFIM{O}
    op_sampler::O
    J::M
    hx::V
    Ns::Int
    energy_shift::Float64
end

const ISING_OP_SIZE = 5
###############################################################################

# TFIM ops:
#  (1,-2,i,0,i) is an off-diagonal site operator h*sigma^x_i
#  (1,1,i,0,i) is a diagonal site operator h
#  (0,0,0,0,0) is the identity operator I - NOT USED IN THE PROJECTOR CASE
#  (2,1,w,i,j) is a diagonal bond operator J(sigma^z_i sigma^z_j)
@inline getoperatorlocality(::Type{<:AbstractIsing}, op::NTuple{ISING_OP_SIZE, Int}) = @inbounds op[1]
@inline getoperatorlocality(H::AbstractIsing, op::NTuple{ISING_OP_SIZE, Int}) = getoperatorlocality(typeof(H), op)
@inline getoperatortype(::Type{<:AbstractIsing}, op::NTuple{ISING_OP_SIZE, Int}) = @inbounds op[2]
@inline getoperatortype(H::AbstractIsing, op::NTuple{ISING_OP_SIZE, Int}) = getoperatortype(typeof(H), op)
@inline getweightindex(::Type{<:AbstractIsing}, op::NTuple{ISING_OP_SIZE, Int}) = @inbounds op[3]
@inline getweightindex(H::AbstractIsing, op::NTuple{ISING_OP_SIZE, Int}) = getweightindex(typeof(H), op)
@inline getsite(::Type{<:AbstractIsing}, op::NTuple{ISING_OP_SIZE, Int}) = @inbounds op[end]
@inline getsite(H::AbstractIsing, op::NTuple{ISING_OP_SIZE, Int}) = getsite(typeof(H), op)
@inline getbondsites(::Type{<:AbstractIsing}, op::NTuple{ISING_OP_SIZE, Int}) = @inbounds (op[end-1], op[end])
@inline getbondsites(H::AbstractIsing, op::NTuple{ISING_OP_SIZE, Int}) = getbondsites(typeof(H), op)

@inline convertoperatortype(H::Type{<:AbstractIsing}, op::NTuple{ISING_OP_SIZE, Int}, new_t::Int) =
    @inbounds (op[1], new_t, getweightindex(H, op) - getoperatortype(H, op) + new_t, op[4], op[5])
@inline convertoperatortype(H::AbstractIsing, op::NTuple{ISING_OP_SIZE, Int}, new_t::Int) = convertoperatortype(typeof(H), op, new_t)

@inline isidentity(H::Type{<:AbstractIsing}, op::NTuple{ISING_OP_SIZE, Int}) = (getoperatorlocality(H, op) == 0)
@inline isidentity(H::AbstractIsing, op::NTuple{ISING_OP_SIZE, Int}) = isidentity(typeof(H), op)
@inline issiteoperator(H::Type{<:AbstractIsing}, op::NTuple{ISING_OP_SIZE, Int}) = (getoperatorlocality(H, op) == 1)
@inline issiteoperator(H::AbstractIsing, op::NTuple{ISING_OP_SIZE, Int}) = issiteoperator(typeof(H), op)
@inline isbondoperator(H::Type{<:AbstractIsing}, op::NTuple{ISING_OP_SIZE, Int}) = (getoperatorlocality(H, op) == 2)
@inline isbondoperator(H::AbstractIsing, op::NTuple{ISING_OP_SIZE, Int}) = isbondoperator(typeof(H), op)
@inline isdiagonal(H::Type{<:AbstractIsing}, op::NTuple{ISING_OP_SIZE, Int}) = (getoperatortype(H, op) != -2)
@inline isdiagonal(H::AbstractIsing, op::NTuple{ISING_OP_SIZE, Int}) = isdiagonal(typeof(H), op)

@inline makeidentity(::Type{<:AbstractIsing})::NTuple{ISING_OP_SIZE, Int} = (0, 0, 0, 0, 0)
@inline makediagonalsiteop(::Type{<:AbstractIsing}, i::Int)::NTuple{ISING_OP_SIZE, Int} = (1, 1, i, 0, i)
@inline makeoffdiagonalsiteop(::Type{<:AbstractIsing}, i::Int)::NTuple{ISING_OP_SIZE, Int} = (1, -2, i, 0, i)
@inline makeidentity(H::AbstractIsing) = makeidentity(typeof(H))
@inline makediagonalsiteop(H::AbstractIsing, i::Int) = makediagonalsiteop(typeof(H), i)
@inline makeoffdiagonalsiteop(H::AbstractIsing, i::Int) = makeoffdiagonalsiteop(typeof(H), i)

@inline getbondtype(::AbstractTFIM, s1::Bool, s2::Bool) = 1

@inline diagonaloperator(::Type{<:AbstractIsing}) = Diagonal([-1, 1])
@inline diagonaloperator(H::AbstractIsing) = diagonaloperator(typeof(H))

@inline getlogweight(H::AbstractIsing, op::NTuple{ISING_OP_SIZE, Int}) = @inbounds getlogweight(H.op_sampler, getweightindex(H, op))

###############################################################################

function make_prob_vector(J::UpperTriangular{T}, hx::AbstractVector{T}) where T
    @assert length(hx) == size(J, 1) == size(J, 2)

    ops = Vector{NTuple{ISING_OP_SIZE, Int}}(undef, 0)
    p = Vector{T}(undef, 0)
    energy_shift = zero(T)

    for i in eachindex(hx)
        if !iszero(hx[i])
            push!(ops, makediagonalsiteop(AbstractTFIM, i))
            push!(p, hx[i])
            energy_shift += hx[i]
        end
    end

    # only take J_ij terms from upper triangle
    for j in axes(J, 2), i in axes(J, 1)
        if i < j  # i != j: we don't want self-interactions
            if J[i, j] != 0
                push!(p, 2*abs(J[i, j]))
                push!(ops, (2, 1, length(p), i, j))
                energy_shift += abs(J[i, j])
            end
        end
    end

    return ops, p, energy_shift
end

function make_uniform_tfim(bond_spins::Vector{NTuple{2,Int}}, Ns::Int, J::T, hx::T) where T
    hx_ = hx * ones(T, Ns)
    J_ = zeros(T, Ns, Ns)
    for (i, j) in bond_spins
        i, j = (i <= j) ? (i, j) : (j, i)
        J_[i, j] = -J
    end

    return UpperTriangular(triu!(J_, 1)), hx_
end

###############################################################################

function TFIM(J::UpperTriangular{Float64}, hx::AbstractVector{Float64})
    @assert length(hx) == size(J, 1) == size(J, 2)

    ops, p, energy_shift = make_prob_vector(J, hx)
    Ns = length(hx)
    op_sampler = OperatorSampler(ops, p)

    return TFIM{typeof(op_sampler), typeof(J), typeof(hx)}(
        op_sampler, J, hx, Ns, energy_shift
    )
end


function TFIM(bond_spin, Ns::Int, Nb::Int, hx::Float64, J::Float64)
    J_, hx_ = make_uniform_tfim(bond_spin, Ns, J, hx)
    return TFIM(J_, hx_)
end

abstract type HXField; end
struct ConstantHX <: HXField; end
struct VaryingHX <: HXField; end

hxfield(::AbstractIsing) = VaryingHX()
hxfield(::NearestNeighbourTFIM) = ConstantHX()

total_hx(::ConstantHX, H::AbstractIsing) = nspins(H) * H.hx
total_hx(::VaryingHX, H::AbstractIsing) = sum(H.hx)
total_hx(H::AbstractIsing) = total_hx(hxfield(H), H)

# returns true if H.J[site1, site2] is negative
Base.@propagate_inbounds isferromagnetic(H::TFIM, (site1, site2)::NTuple{2, Int}) = signbit(H.J[site1, site2])
haslongitudinalfield(::AbstractTFIM) = false


###############################################################################


function energy(::Type{<:BinaryGroundState}, H::AbstractIsing, ns::Vector{<:Real}; resampler::Function=jackknife)
    hx = total_hx(H)

    if !iszero(hx)
        E = -hx * resampler(inv, ns)
    else
        E = measurement(zero(H.energy_shift))
    end

    return H.energy_shift + E
end

function energy(::Type{<:BinaryGroundState}, H::AbstractIsing, ns_mean::Real)
    hx = total_hx(H)

    if !iszero(hx)
        E = -hx / ns_mean
    else
        E = zero(H.energy_shift)
    end

    return H.energy_shift + E
end
