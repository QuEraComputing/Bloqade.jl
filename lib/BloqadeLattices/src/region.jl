abstract type AbstractRegion{D} end
####### AbstractRegion ########
# defines a finite region in R^D which will be used to generate bounded lattice_sites
# The region must have the following functions
# 

struct Parallelepiped{D, T} <: AbstractRegion{D}
    # abstract matrix dimensions already guaranteed to be 2D
    vecs::Matrix{T}
    vecs_inv::Matrix{T}
end


function Parallelepiped(vecs)
    D = size(vecs,1)
    vecs_inv = inv(vecs)
    # return Parallelpiped
    return Parallelepiped{D,eltype(vecs_inv)}(vecs,vecs_inv)
end

# Handle the 1D case
function Parallelepiped(vecs::T) where {T<:Real}
    vecs = fill(vecs,1,1)
    vecs_inv = inv(vecs)
    # return Parallelpiped
    return Parallelepiped{1,eltype(vecs_inv)}(vecs,vecs_inv)
end


Base.broadcastable(x::AbstractRegion) = Ref(x)

# check if a point is in the region
in_range(x) = 0 ≤ x < 1 ? true : false
Base.in(x,region::Parallelepiped{D,T}) where {D,T} = all(in_range.(region.vecs_inv * [x...,]))

# Enforce periodic boundary conditions by having points that fall outside of the 
# parallelogram map to ones on the inside
Base.mod(x,region::Parallelepiped{D,T}) where {D,T} = typeof(x)(region.vecs * mod.((region.vecs_inv * [x...,]),1.0))
Base.mod(x,region::Parallelepiped{1,T}) where {T} =  mod.(x,only(region.vecs))

distance(x, y) = sqrt(mapreduce(x -> x^2, +, x .- y))

function distance(region::Parallelepiped{D,T},x,y) where {D,T}
    x = mod([x...,],region)
    y = mod([y...,],region)
    dist = distance(x,y)

    for a in 1:D
        shift = region.vecs[:,a]
        for b in a:D
            shift .+= (b > a ? region.vecs[:,b] : 0)
            shift_y = y .+ shift
            shift_x = x .+ shift
            dist = min(
                distance(x,shift_y),
                distance(shift_x,y),
                dist
            )
        end
    end

    return dist
end

