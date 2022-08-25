abstract type Lattice end
abstract type OneDLattice <: Lattice end
abstract type BravaisLattice <: Lattice end
abstract type PolyLattice <: Lattice end

# 2D and 1D lattices
# 1D: even-spaced chains only (for now)
# 2D: Simple lattice types are derived from modifying a monoclinic lattice's parameters
# Special lattices like Ruby, Kagome, are custom-made

# SimpleLattice: A lattice with a one-site basis (e.g. rectangular, 1D chain)
# PolyLattice: A lattice with a multi-site basis, but is common (e.g. Kagome)
# TODO: 1D chains with multi-site basis, other common lattices


# 1D chain

struct Chain <: OneDLattice
    n1::Int
    a1::Float64
    distance_matrix::Array{Float64, 2}
end

function Chain(nX::Int, aX::Float64, PBC::Bool; trunc::Float64 = Inf)
    distance_matrix = zeros(Float64, nX, nX)

    for i in 1:(nX-1)
        x1 = rem(i, nX) > 0 ? rem(i, nX) : nX

        for j in (i+1):nX
            x2 = rem(j, nX) > 0 ? rem(j, nX) : nX

            if PBC
                dx = abs(x1-x2) > 0.5*nX ? 0.5*nX - rem(abs(x1-x2), 0.5*nX) : abs(x1-x2)
            else
                dx = abs(x1-x2)
            end

            dx *= aX
            distance_matrix[i, j] = dx <= trunc ? dx : 0.0
        end
    end
    return Chain(nX, aX, distance_matrix)
end


# Monoclinic lattice

struct Monoclinic <: BravaisLattice
    n1::Int
    n2::Int
    a::Float64
    b::Float64
    θ::Float64
    PBC::Tuple{Bool, Bool}
    distance_matrix::Array{Float64, 2}
end

function Monoclinic(n1::Int, n2::Int, a::Float64, b::Float64, θ::Float64, PBC::Tuple{Bool, Bool}; trunc::Float64 = Inf)
    # primitive lattice translation vectors:
    a1 = [a, 0.0]
    a2 = b * [cos(θ), sin(θ)]
    distance_matrix = create_simple_distance_matrix(n1, n2, a1, a2, PBC; trunc=trunc)
    return Monoclinic(n1, n2, a1, a2, PBC, distance_matrix)
end
Monoclinic(n1::Int, n2::Int, a::Float64, b::Float64, θ::Float64, PBC::Bool; trunc::Float64 = Inf) = Monoclinic(n1, n2, a, b, θ, (PBC, PBC); trunc = trunc)


# Rectangular

struct Rectangle <: BravaisLattice
    n1::Int
    n2::Int
    a::Float64
    b::Float64
    PBC::Tuple{Bool, Bool}
    distance_matrix::Array{Float64, 2}
end

function Rectangle(n1::Int, n2::Int, a::Float64, b::Float64, PBC::Tuple{Bool, Bool}; trunc::Float64 = Inf)
    a1 = [a, 0.0]
    a2 = [0.0, b]
    distance_matrix = create_simple_distance_matrix(n1, n2, a1, a2, PBC; trunc=trunc)
    return Rectangle(n1, n2, a, b, PBC, distance_matrix)
end
Rectangle(n1::Int, n2::Int, a::Float64, b::Float64, PBC::Bool; trunc::Float64 = Inf) = Rectangle(n1, n2, a, b, (PBC, PBC); trunc = trunc)


# Square

struct Square <: BravaisLattice
    n1::Int
    n2::Int
    a::Float64
    PBC::Tuple{Bool, Bool}
    distance_matrix::Array{Float64, 2}
end

function Square(n1::Int, n2::Int, a::Float64, PBC::Tuple{Bool, Bool}; trunc::Float64 = Inf)
    a1 = [a, 0.0]
    a2 = [0.0, a]
    distance_matrix = create_simple_distance_matrix(n1, n2, a1, a2, PBC; trunc=trunc)
    return Square(n1, n2, a, PBC, distance_matrix)
end
Square(n1::Int, n2::Int, a::Float64, PBC::Bool; trunc = Inf) = Square(n1, n2, a, a, (PBC, PBC))


# Honeycomb

struct Honeycomb <: BravaisLattice
    n1::Int
    n2::Int
    a::Float64
    b::Float64
    PBC::Tuple{Bool, Bool}
    distance_matrix::Array{Float64, 2}
end

function Honeycomb(n1::Int, n2::Int, a::Float64, PBC::Tuple{Bool, Bool}; trunc::Float64 = Inf)
    a1 = [a, 0.0]
    a2 = a .* [cos(2*π/3.), sin(2*π/3.)]
    distance_matrix = create_simple_distance_matrix(n1, n2, a1, a2, PBC; trunc=trunc)
    return Rectangle(n1, n2, a, b, PBC, distance_matrix)
end
Honeycomb(n1::Int, n2::Int, a::Float64, PBC::Bool; trunc::Float64 = Inf) = Honeycomb(n1, n2, a, (PBC, PBC); trunc = trunc)


###############################################################################


# Kagome lattice

struct Kagome <: PolyLattice
    # number of repititions in directions of a1 and a2
    n1::Int # a1
    n2::Int # a2
    # parameter that defines the equilateral triangle side length
    t::Float64
    # translation vectors
    a1::Array{Float64,1}
    a2::Array{Float64,1}
    # coordinates of sites inside a unit cell
    r::Array{Array{Float64,1},1}

    PBC::Tuple{Bool,Bool}
    distance_matrix::Array{Float64, 2}
end

function Kagome(t::Float64, n1::Int, n2::Int, PBC::Tuple{Bool, Bool}; trunc::Float64 = Inf)
    a = 2*t
    a1 = [a, 0.]
    a2 = [a*0.5, a*sqrt(3)*0.5]

    # coordinates of each site in the unit cell
    r1 = [0., 0.]
    r2 = 0.5 * a2
    r3 = 0.5 * a1
    r = [r1, r2, r3]

    distance_matrix = create_distance_matrix(n1, n2, a1, a2, r, PBC; trunc=trunc)

    return Kagome(n1, n2, t, a1, a2, r, PBC, distance_matrix)
end
Kagome(t::Float64, n1::Int, n2::Int, PBC::Bool; trunc::Float64 = Inf) = Kagome(t, n1, n2, (PBC, PBC); trunc = trunc)


# Ruby lattice

struct Ruby <: PolyLattice
    # number of repititions in directions of a1 and a2
    n1::Int # a1
    n2::Int # a2
    # parameter that defines spacing
    ρ::Float64
    # translation vectors
    a1::Array{Float64,1}
    a2::Array{Float64,1}
    # coordinates of sites inside a unit cell
    r::Array{Array{Float64,1},1}

    PBC::Tuple{Bool,Bool}
    distance_matrix::Array{Float64, 2}
end

function Ruby(ρ::Float64, n1::Int, n2::Int, PBC::Tuple{Bool, Bool}; trunc::Float64 = Inf)
    a = 4*ρ/sqrt(3)
    a1 = [a, 0.]
    a2 = [a*0.5, a*sqrt(3)*0.5]

    # coordinates of each site in the unit cell
    r1 = [0., 0.]
    r2 = 0.75 * a2
    r3 = 0.25 * (a1 + a2)
    r4 = 0.5 * a1
    r5 = 0.25*a1 + 0.75*a2
    r6 = 0.5*a1 + 0.25*a2
    r = [r1, r2, r3, r4, r5, r6]

    distance_matrix = create_distance_matrix(n1, n2, a1, a2, r, PBC; trunc=trunc)

    return Ruby(n1, n2, ρ, a1, a2, r, PBC, distance_matrix)
end
Ruby(ρ::Float64, n1::Int, n2::Int, PBC::Bool; trunc::Float64 = Inf) = Ruby(ρ, n1, n2, (PBC, PBC); trunc = trunc)


# Custom lattice. Use if required lattice is not specified above

struct Custom <: PolyLattice
    # number of repititions in directions of a1 and a2
    n1::Int # a1
    n2::Int # a2
    # translation vectors
    a1::Array{Float64,1}
    a2::Array{Float64,1}
    # coordinates of sites inside a unit cell
    r::Array{Array{Float64,1},1}

    # Periodic boundaries in directions of a1, a2
    # these are conventional PBCs
    PBC::Tuple{Bool,Bool}
    distance_matrix::Array{Float64, 2}
end

function Custom(n1::Int, n2::Int, a1::Vector{Float64}, a2::Vector{Float64}, r::Array{Array{Float64,1},1}, PBC::Tuple{Bool, Bool}; trunc::Float64 = Inf)
    distance_matrix = create_distance_matrix(n1, n2, a1, a2, r, PBC; trunc=trunc)
    return Custom(n1, n2, a1, a2, r, PBC, distance_matrix)
end

nspins(lattice::OneDLattice) = lattice.n1
nspins(lattice::BravaisLattice) = lattice.n1 * lattice.n2
nspins(lattice::PolyLattice) = size(lattice.r)[1] * lattice.n1 * lattice.n2


###############################################################################

#=
function circleShape(h, k, r)
    θ = LinRange(0, 2*π, 500)
    h .+ r*sin.(θ), k .+ r*cos.(θ)
end
=#

# Distance matrix functions

function create_simple_distance_matrix(n1::Int, n2::Int, a1::Vector{Float64}, a2::Vector{Float64}, PBC::Tuple{Bool, Bool}; trunc::Float64 = Inf)
    PBC1, PBC2 = PBC

    a2_length = sqrt(a2[1]^2 + a2[2]^2) # length of a2 vector
    θ = acos(a2[1] / a2_length) # a2 angle from horizontal

    # total number of sites
    N = n1*n2
    dij = zeros(Float64, N, N)

    # (xpi, ypi) are lattice index coordinates
    # need to keep track of these for PBC enforcement
    # subtract one from xpi and ypi to center everything at (0,0)
    for i in 1:(N-1)
        xp1 = rem(i, n1) > 0 ? rem(i, n1) - 1 : n1 - 1
        yp1 = rem(i, n1) > 0 ? div(i, n1) : div(i, n1) - 1

        x1 = xp1*a1[1] + yp1*a2_length*cos(θ)
        y1 = yp1*a2_length*sin(θ)

        for j in (i+1):N
            xp2 = rem(j, n1) > 0 ? rem(j, n1) - 1 : n1 - 1
            yp2 = rem(j, n1) > 0 ? div(j, n1) : div(j, n1) - 1

            x2 = xp2*a1[1] + yp2*a2_length*cos(θ)
            y2 = yp2*a2_length*sin(θ)

            # check x and y directions to enforce PBCs
            if PBC1 & PBC2
                # calculate minimum distance

                # non-periodic
                d_np = sqrt( (y1 - y2)^2 + (x1 - x2)^2 )

                # periodic in x only
                x2 -= a1[1]*n1*cos(θ)
                d_px = sqrt( (y1 - y2)^2 + (x1 - x2)^2 )

                # periodic in x and y
                # x2 has already been taken care of
                y2 -= a2_length*n2*sin(θ)
                d_pxy = sqrt( (y1 - y2)^2 + (x1 - x2)^2 )

                # periodic in y only
                # undo the x2 periodicity
                x2 += a1[1]*n1*cos(θ)
                d_py = sqrt( (y1 - y2)^2 + (x1 - x2)^2 )

                d = min(d_py, d_px, d_pxy, d_np)

            elseif PBC1 & !PBC2
                if abs(x1 - x2) > 0.5*n1*a1[1]
                    x2 -= a1[1]*n1*cos(θ)
                end

                d = sqrt( (y1 - y2)^2 + (x1 - x2)^2 )

            elseif PBC2 & !PBC1
                yi = x1*cos(θ) + y1*sin(θ)
                yj = x2*cos(θ) + y2*sin(θ)

                if abs(yi - yj) > 0.5*n2*a2_length
                    x2 -= a2_length*n2*cos(θ)
                    y2 -= a2_length*n2*sin(θ)
                end

                d = sqrt( (y1 - y2)^2 + (x1 - x2)^2 )

            else # no PBCs anywhere
                d = sqrt( (y1 - y2)^2 + (x1 - x2)^2 )
            end

            dij[i, j] = d <= trunc ? d : 0.0
        end
    end

    return dij
end


function create_distance_matrix(n1::Int, n2::Int, a1::Vector{Float64}, a2::Vector{Float64}, r::Array{Array{Float64,1},1}, PBC::Tuple{Bool, Bool}; trunc::Float64 = Inf, plotname="lattice", plotting=false)
    PBC1, PBC2 = PBC

    a2_length = sqrt(a2[1]^2 + a2[2]^2) # length of a2 vector
    θ = acos(a2[1] / a2_length) # a2 angle from horizontal

    N = n1*n2*length(r) # total number of sites
    dij = zeros(Float64, N, N) # distance matrix
    num_cells = n1*n2 # number of repitions of the unit cell

    a2_length = sqrt(a2[1]^2 + a2[2]^2) # length of a2 vector
    θ = acos(a2[1] / a2_length) # a2 angle from horizontal

    if plotting
        xs = []
        ys = []
    end
    
    # NOT ENDING ON num_cells-1 BECAUSE WE NEED INTER-CELL BONDS
    for i in 1:num_cells
        cellp1_x = rem(i, n1) > 0 ? rem(i, n1) - 1 : n1 - 1
        cellp1_y = rem(i, n1) > 0 ? div(i, n1) : div(i, n1) - 1

        cell1_x = cellp1_x*a1[1] + cellp1_y*a2_length*cos(θ)
        cell1_y = cellp1_y*a2_length*sin(θ)
        cell1   = [cell1_x, cell1_y]

        # NOT STARTING FROM i+1 BECAUSE WE NEED INTER-CELL BONDS
        for j in i:num_cells
            cellp2_x = rem(j, n1) > 0 ? rem(j, n1) - 1 : n1 - 1
            cellp2_y = rem(j, n1) > 0 ? div(j, n1) : div(j, n1) - 1

            cell2_x = cellp2_x*a1[1] + cellp2_y*a2_length*cos(θ)
            cell2_y = cellp2_y*a2_length*sin(θ)
            cell2   = [cell2_x, cell2_y]

            for site_i in 1:length(r)
                site_num_i = site_i + length(r)*(i-1)
                ri = r[site_i] + cell1
                
                if plotting
                    push!(xs, ri[1])
                    push!(ys, ri[2])
                end
                
                for site_j in 1:length(r)
                    site_num_j = site_j + length(r)*(j-1)

                    # ensure that there are no diagonal entries
                    if site_num_i == site_num_j
                        continue
                    end

                    rj = r[site_j] + cell2

                    #=
                    if plotting
                        push!(xs, rj[1])
                        push!(ys, rj[2])
                    end
                    =#

                    Δ = ri - rj
                    Δx, Δy = Δ[1], Δ[2]
                    d_np = sqrt(Δx^2 + Δy^2)

                    # checks for PBCs
                    if PBC1 & PBC2
                        # periodic in both lattice directions
                        # essentially, find smallest distance

                        # distances for periodicity in a1
                        rj .+= a1 * n1
                        Δ = ri - rj
                        Δx, Δy = Δ[1], Δ[2]
                        d_1 = sqrt(Δx^2 + Δy^2)

                        rj .-= 2 * a1 * n1
                        Δ = ri - rj
                        Δx, Δy = Δ[1], Δ[2]
                        d_2 = sqrt(Δx^2 + Δy^2)

                        # undo change
                        rj .+= a1 * n1

                        # distances for periodicity in a2
                        rj .+= a2 * n1
                        Δ = ri - rj
                        Δx, Δy = Δ[1], Δ[2]
                        d_3 = sqrt(Δx^2 + Δy^2)

                        rj .-= 2 * a2 * n1
                        Δ = ri - rj
                        Δx, Δy = Δ[1], Δ[2]
                        d_4 = sqrt(Δx^2 + Δy^2)

                        # undo change
                        rj .+= a2 * n1

                        # take the minimum distance
                        d = min(d_np, d_1, d_2, d_3, d_4)

                    elseif PBC1 & !PBC2
                        # distances for periodicity in a1
                        rj .+= a1 * n1
                        Δ = ri - rj
                        Δx, Δy = Δ[1], Δ[2]
                        d_1 = sqrt(Δx^2 + Δy^2)

                        rj .-= 2 * a1 * n1
                        Δ = ri - rj
                        Δx, Δy = Δ[1], Δ[2]
                        d_2 = sqrt(Δx^2 + Δy^2)

                        # undo change
                        rj .+= a1 * n1

                        d = min(d_np, d_1, d_2)

                    elseif PBC2 & !PBC1
                        # distances for periodicity in a2
                        rj .+= a2 * n1
                        Δ = ri - rj
                        Δx, Δy = Δ[1], Δ[2]
                        d_3 = sqrt(Δx^2 + Δy^2)

                        rj .-= 2 * a2 * n1
                        Δ = ri - rj
                        Δx, Δy = Δ[1], Δ[2]
                        d_4 = sqrt(Δx^2 + Δy^2)

                        # undo change
                        rj .+= a2 * n1

                        # take the minimum distance
                        d = min(d_np, d_3, d_4)

                    else
                        d = d_np
                    end

                    dij[site_num_i, site_num_j] = d <= trunc ? d : 0.0

                end
            end
        end
    end

    if plotting
        plot(xs, ys, seriestype = :scatter)
        savefig(plotname * ".png")
        points = hcat(xs, ys)
        writedlm("coords.CSV", points, ',')
    end

    return dij
end
create_simple_distance_matrix(n1::Int, n2::Int, a1::Vector{Float64}, a2::Vector{Float64}, PBC::Bool; trunc::Float64 = Inf) = create_simple_distance_matrix(n1, n2, a1, a2, (PBC, PBC); trunc = trunc)
create_distance_matrix(n1::Int, n2::Int, a1::Vector{Float64}, a2::Vector{Float64}, r::Array{Array{Float64,1},1}, PBC::Bool; trunc::Float64 = Inf) = create_distance_matrix(n1, n2, a1, a2, r, (PBC, PBC); trunc = trunc)