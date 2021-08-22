

using LuxurySparse
using Serialization
atoms = square_lattice(40, 0.8)
graph = unit_disk_graph(atoms, 1.5)
space = blockade_subspace(graph)

# write_atoms("test.atoms", atoms)
# atoms = read_atoms("test.atoms")
# space = deserialize("subspace.dat")


h = RydInteract(atoms) + XTerm(length(atoms), 1.0) - NTerm(length(atoms), 1.2)

# SparseMatrixCSC{ComplexF32}(RydInteract(atoms), space)
SparseMatrixCSC{ComplexF32}(h, space)
# SparseMatrixCSC{ComplexF32}(NTerm(length(atoms), 1.2), space)


H = SparseMatrixCOO{ComplexF32}(undef, length(space), length(space));
to_matrix!(H, XTerm(length(atoms), 1.0), space);

maximum(H.is)
maximum(H.js)
serialize("test_mat.dat", H)
S = SparseMatrixCSC(H)
