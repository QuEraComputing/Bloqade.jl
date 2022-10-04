abstract type AbstractSuperCell{D} end

struct Parallelepiped{D, T} <: AbstractSuperCell{D}
    # abstract matrix dimensions already guaranteed to be 2D
    bounds::Matrix{T}
    bounds_inv::Matrix{T}
end


function Parallelepiped(bounds)
    D = size(bounds,1)
    bounds_inv = inv(bounds)
    # return Parallelpiped
    return Parallelepiped{D,eltype(bounds_inv)}(bounds,bounds_inv)
end

# Handle the 1D case
function Parallelepiped(bounds::T) where {T<:Real}
    bounds = fill(bounds,1,1)
    bounds_inv = inv(bounds)
    # return Parallelpiped
    return Parallelepiped{1,eltype(bounds_inv)}(bounds,bounds_inv)
end


Base.broadcastable(x::AbstractSuperCell) = Ref(x)

# check if a point is in the cell
in_range(x) = 0 ≤ x < 1 ? true : false
within_cell(cell::Parallelepiped{D,T},x) where {D,T} = all(in_range.(cell.bounds_inv * [x...,]))

# Enforce periodic boundary conditions by having points that fall outside of the 
# parallelogram map to ones on the inside
wrap_around(cell::Parallelepiped{D,T},x) where {D,T} = typeof(x)(cell.bounds * mod.((cell.bounds_inv * [x...,]),1.0))
wrap_around(cell::Parallelepiped{1,T},x) where {T} =  mod.(x,only(cell.bounds))

distance(x, y) = sqrt(mapreduce(x -> x^2, +, x .- y))

function distance(cell::Parallelepiped{D,T},x,y) where {D,T}
    x = wrap_around(cell,[x...,])
    y = wrap_around(cell,[y...,])
    dist = distance(x,y)

    for a in 1:D
        shift = cell.bounds[:,a]
        for b in a:D
            shift = shift .+ (b > a ? cell.bounds[:,b] : 0)
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

