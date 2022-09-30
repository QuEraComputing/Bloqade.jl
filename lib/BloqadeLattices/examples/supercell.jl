using BloqadeLattices
using Base.Iterators

function wrap_around(cell::Parallelepiped{D,T},x) where {D,T}
    return reshape(cell.bounds * ((cell.bounds_inv * [x...,] .+ one(T)) .% 1) , size(x)...)
end


function distance(cell::Parallelepiped{D,T},x,y) where {D,T}
    x = wrap_around(cell,[x...,])
    y = wrap_around(cell,[y...,])
    dist = BloqadeLattices.distance(x,y)

    for a in 1:D
        shift = cell.bounds[:,a]
        for b in a:D
            shift .+= (b > a ? cell.bounds[:,b] : 0)
            shift_y = y .+ shift
            shift_x = x .+ shift
            dist = min(
                BloqadeLattices.distance(x,shift_y),
                BloqadeLattices.distance(shift_x,y),
                dist
            )
        end
    end

    return dist
end

begin
    bounds = zeros((2,2))
    bounds[:,1] .= (1,3)
    bounds[:,2] .= (3,0)

    pp = Parallelepiped(bounds)

    x = (1,1)
    y = (3,2)

    println(distance(pp,y,x)," ",distance(pp,x,y))



end

