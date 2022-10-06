###############################################################################

# TFIM ops:
#  (1,-2,i,0,i) is an off-diagonal site operator h*sigma^x_i
#  (1,1,i,0,i) is a diagonal site operator h
#  (0,0,0,0,0) is the identity operator I - NOT USED IN THE PROJECTOR CASE
#  (2,1,w,i,j) is a diagonal bond operator J(sigma^z_i sigma^z_j)


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

# from previous Rydberg.jl:

@inline diagonaloperator(::Type{<:AbstractRydberg}) = Diagonal([0, 1])
@inline diagonaloperator(H::AbstractRydberg) = diagonaloperator(typeof(H))
