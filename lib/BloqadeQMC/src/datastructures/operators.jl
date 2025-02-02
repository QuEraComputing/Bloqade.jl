###############################################################################

# Local Rydberg ops:
#  (1,-2,i,0,i) is an off-diagonal site operator h*sigma^x_i
#  (1,1,i,0,i) is a diagonal site operator h
#  (0,0,0,0,0) is the identity operator I - NOT USED IN THE PROJECTOR CASE
#  (2,1,w,i,j) is a diagonal bond operator J(sigma^z_i sigma^z_j)

#   The weight index encodes the spin configuration between which the operator is sandwiched.  It encodes this by subtracting 1, then converting to binary. 
#    w = 1 -> 00 -> down-down
#    w = 2 -> 01 -> down-up
#    w = 3 -> 10 -> up-down
#    w = 4 -> 11 -> up-up
#    See getbondtype() for implementation. (That function is called in the multibranch and line kernel functions.)

###############################################################################

OP_SIZE = 5     # number of indices used to specify a local operator

@inline getbondtype(::AbstractRydberg, s1::Bool, s2::Bool) = (s1<<1 | s2) + 1
@inline spin_config(::AbstractRydberg, t::Int)::NTuple{2,Int} = divrem(t - 1, 2)
@inline spin_config(H::AbstractRydberg, op::NTuple{OP_SIZE, Int}) = @inbounds spin_config(H, getoperatortype(H, op))

###############################################################################

@inline getoperatorlocality(::Type{<:AbstractRydberg}, op::NTuple{OP_SIZE, Int}) = @inbounds op[1]
@inline getoperatorlocality(H::AbstractRydberg, op::NTuple{OP_SIZE, Int}) = getoperatorlocality(typeof(H), op)
@inline getoperatortype(::Type{<:AbstractRydberg}, op::NTuple{OP_SIZE, Int}) = @inbounds op[2]
@inline getoperatortype(H::AbstractRydberg, op::NTuple{OP_SIZE, Int}) = getoperatortype(typeof(H), op)
@inline getweightindex(::Type{<:AbstractRydberg}, op::NTuple{OP_SIZE, Int}) = @inbounds op[3]
@inline getweightindex(H::AbstractRydberg, op::NTuple{OP_SIZE, Int}) = getweightindex(typeof(H), op)
@inline getsite(::Type{<:AbstractRydberg}, op::NTuple{OP_SIZE, Int}) = @inbounds op[end]
@inline getsite(H::AbstractRydberg, op::NTuple{OP_SIZE, Int}) = getsite(typeof(H), op)
@inline getbondsites(::Type{<:AbstractRydberg}, op::NTuple{OP_SIZE, Int}) = @inbounds (op[end-1], op[end])
@inline getbondsites(H::AbstractRydberg, op::NTuple{OP_SIZE, Int}) = getbondsites(typeof(H), op)

@inline convertoperatortype(H::Type{<:AbstractRydberg}, op::NTuple{OP_SIZE, Int}, new_t::Int) =
    @inbounds (op[1], new_t, getweightindex(H, op) - getoperatortype(H, op) + new_t, op[4], op[5])
@inline convertoperatortype(H::AbstractRydberg, op::NTuple{OP_SIZE, Int}, new_t::Int) = convertoperatortype(typeof(H), op, new_t)

@inline isidentity(H::Type{<:AbstractRydberg}, op::NTuple{OP_SIZE, Int}) = (getoperatorlocality(H, op) == 0)
@inline isidentity(H::AbstractRydberg, op::NTuple{OP_SIZE, Int}) = isidentity(typeof(H), op)
@inline issiteoperator(H::Type{<:AbstractRydberg}, op::NTuple{OP_SIZE, Int}) = (getoperatorlocality(H, op) == 1)
@inline issiteoperator(H::AbstractRydberg, op::NTuple{OP_SIZE, Int}) = issiteoperator(typeof(H), op)
@inline isbondoperator(H::Type{<:AbstractRydberg}, op::NTuple{OP_SIZE, Int}) = (getoperatorlocality(H, op) == 2)
@inline isbondoperator(H::AbstractRydberg, op::NTuple{OP_SIZE, Int}) = isbondoperator(typeof(H), op)
@inline isdiagonal(H::Type{<:AbstractRydberg}, op::NTuple{OP_SIZE, Int}) = (getoperatortype(H, op) != -2)
@inline isdiagonal(H::AbstractRydberg, op::NTuple{OP_SIZE, Int}) = isdiagonal(typeof(H), op)

@inline makeidentity(::Type{<:AbstractRydberg})::NTuple{OP_SIZE, Int} = (0, 0, 0, 0, 0)
@inline makediagonalsiteop(::Type{<:AbstractRydberg}, i::Int)::NTuple{OP_SIZE, Int} = (1, 1, i, 0, i)
@inline makeoffdiagonalsiteop(::Type{<:AbstractRydberg}, i::Int)::NTuple{OP_SIZE, Int} = (1, -2, i, 0, i)
@inline makeidentity(H::AbstractRydberg) = makeidentity(typeof(H))
@inline makediagonalsiteop(H::AbstractRydberg, i::Int) = makediagonalsiteop(typeof(H), i)
@inline makeoffdiagonalsiteop(H::AbstractRydberg, i::Int) = makeoffdiagonalsiteop(typeof(H), i)

@inline getbondtype(::AbstractTFIM, s1::Bool, s2::Bool) = 1

@inline diagonaloperator(::Type{<:AbstractRydberg}) = Diagonal([-1, 1])
@inline diagonaloperator(H::AbstractRydberg) = diagonaloperator(typeof(H))

@inline getlogweight(H::AbstractRydberg, op::NTuple{OP_SIZE, Int}) = @inbounds getlogweight(H.op_sampler, getweightindex(H, op))

###############################################################################

# from previous Rydberg.jl:

@inline diagonaloperator(::Type{<:AbstractRydberg}) = Diagonal([0, 1])
@inline diagonaloperator(H::AbstractRydberg) = diagonaloperator(typeof(H))
