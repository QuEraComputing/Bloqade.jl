Base.@kwdef struct LatticeRow{T <: Real}
    xs::Vector{T} # all the x coordinates for atoms within a row unit cell.
    x_unit_cell_distance::T # the distance between two nearest neighbor x unit cells.
    nrepetitions::Int # the number of repetitions of the unit cell.
end

Base.@kwdef struct FangliLattice{T <: Real}
    rows::Vector{LatticeRow{T}}
    ys::Vector{T} # y is a vector that contains the y coordinates for each row.
    y_unit_cell_distance::T # unit
    nrepetitions::Int # number of repetitions
    x_shifts::Vector{T} # shift of the big unit cell on x axis
end

build_row(row::LatticeRow) = build_row(row.xs, row.x_unit_cell_distance, row.nrepetitions)

# build single row along the x direction
# output: vector containing all x coordinates for atoms in the row 
function build_row(xs::Vector{T}, x_unit_cell_distance::T, nrepetitions::Int) where {T <: Real}
    x_coords = Vector{T}(undef, length(xs) * nrepetitions)

    for i=1:nrepetitions
        for j in 1:length(xs)
            x_coords[(i - 1) * length(xs) + j] = xs[j] + (i - 1) * x_unit_cell_distance
        end
    end

    return x_coords
    
end

# combine several single rows into super-row and build entire 
# lattice by repeating y unit cells 
function site_positions(lattice::FangliLattice{T}) where {T <: Real}
    positions = RydAtom{2, T}[]

    unit_cell_xs = Vector{T}[]
    for row in lattice.rows
        xs = build_row(row)
        push!(unit_cell_xs, xs)
    end
    
    ys = lattice.ys
    for idx in 1:lattice.nrepetitions
        if idx > 1
            ys = ys .+ lattice.y_unit_cell_distance
        end
        shift = lattice.x_shifts[idx]

        for (row_idx, row_xs) in enumerate(unit_cell_xs)
            for col_idx in 1:length(row_xs)
                push!(positions, RydAtom(row_xs[col_idx] + shift, ys[row_idx]))
            end
        end
    end
    return positions
end
