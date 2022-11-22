using BloqadeGates: two_level_indices

@test two_level_indices(1) == [1, 2]
@test two_level_indices(2) == [1, 2, 4, 5]
@test two_level_indices(3) == [1, 2, 4, 5, 10, 11, 13, 14]