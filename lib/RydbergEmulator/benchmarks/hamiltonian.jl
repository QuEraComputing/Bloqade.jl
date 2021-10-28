using SparseArrays
using RydbergEmulator
using BenchmarkTools

atoms = square_lattice(20, 0.8)
h = rydberg_h(atoms, 0.1, nothing, 0.5)
@btime SparseMatrixCSC(h)

lattice = [0 0 1 0 1 1 1;
            1 1 0 1 1 1 0;
            0 1 1 1 1 1 0;
            1 1 1 1 0 1 1;
            0 0 0 0 0 0 0;
            1 1 1 1 0 0 0;
            1 1 1 1 1 0 1]

atoms = [RydAtom(i, j) for i in axes(lattice, 1), j in axes(lattice, 2) if lattice[i, j] == 1]
space = blockade_subspace(atoms, 1.5)
h = rydberg_h(atoms, 0.1, nothing, 0.5)
@btime SparseMatrixCSC(h, space)
