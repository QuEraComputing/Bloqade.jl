using Test
using EaRyd
using EaRyd: build_row, site_positions

right_lattice = [
    FangliLattice(
        rows = [
            LatticeRow([1, 5], 12, 2),
            LatticeRow([3], 12, 2),
        ],
        ys = [1, 3],
        y_unit_cell_distance = 8,
        nrepetitions = 2,
        x_shifts = zeros(Int, 2)
    ),
    FangliLattice(
        rows = [
            LatticeRow([1, 5], 12, 2),
            LatticeRow([3], 12, 2),
        ],
        ys = [1, 3],
        y_unit_cell_distance = 8,
        nrepetitions = 2,
        x_shifts = [-1, 0],
    ),
    FangliLattice(
        rows = [
            LatticeRow([1, 5], 12, 2),
            LatticeRow([3], 12, 2),
        ],
        ys = [1, 5],
        y_unit_cell_distance = 12,
        nrepetitions = 3,
        x_shifts = [-1, 0, 0],
    ),

    FangliLattice(
        rows = [
            LatticeRow([1, 5], 12, 2),
            LatticeRow([3], 3, 2),
        ],
        ys = [1, 5],
        y_unit_cell_distance = 12,
        nrepetitions = 3,
        x_shifts = [-1, 4, 0],
    ),
]

wrong_lattice = [
    FangliLattice(
        rows = [
            LatticeRow([1, 1], 12, 2),
            LatticeRow([3], 3, 2),
        ],
        ys = [1, 5],
        y_unit_cell_distance = 12,
        nrepetitions = 3,
        x_shifts = [-1, 4],
    ),

    FangliLattice(
        rows = [
            LatticeRow([1, 4], 12, 2),
            LatticeRow([3], 3, 2),
        ],
        ys = [3, 3],
        y_unit_cell_distance = 12,
        nrepetitions = 3,
        x_shifts = [-1, 4],
    ),
     
    FangliLattice(
        rows = [
            LatticeRow([1, 4], 12, 2),
            LatticeRow([3], 3, 2),
        ],
        ys = [3, 5, 8],
        y_unit_cell_distance = 12,
        nrepetitions = 3,
        x_shifts = [-1, 4],
    ),
]


@test build_row(right_lattice[1].rows[1]) == [1, 5, 13, 17]
@test build_row(right_lattice[1].rows[2]) == [3, 15]

xs = [1, 5, 13, 17, 3, 15, 1, 5, 13, 17, 3, 15]
ys = [1, 1, 1, 1, 3, 3, 9, 9, 9, 9, 11, 11]

for (idx, pos) in enumerate(site_positions(right_lattice[1]))
    @test xs[idx] == pos[1]
    @test ys[idx] == pos[2]
end

site_positions(right_lattice[2])
site_positions(right_lattice[3])
site_positions(right_lattice[4])