using Test
using BloqadeExpr

positions = rand(2, 10)
@test all(matrix_to_positions(positions)[1] .== positions[:, 1])
