using Test
using Configurations
using OrderedCollections
using EaRydKrylovEvolution

atoms = square_lattice(5, 0.3)
write_atoms("test.atoms", atoms)
@test read_atoms("test.atoms") == atoms

h1 = RydInteract(atoms, 1.0) + XTerm(5, 0.1, 0.2) + ZTerm(5, 0.1)
@test from_toml(Hamiltonian, "hamiltonians/h1.toml") == h1

h2 = XTerm([0.1, 0.2, 0.3, 0.4, 0.5], [0.2, 0.3, 0.4, 0.5, 0.6]) +
    ZTerm([0.1, 0.2, 0.3, 0.4, 0.5]) +
    RydInteract([RydAtom(1, 2) for _ in 1:5], 1.0)

@test from_toml(Hamiltonian, "hamiltonians/h2.toml") == h2

h3 = XTerm([0.1, 0.2, 0.3, 0.4, 0.5], 0.1) + ZTerm(5, 0.1) + RydInteract(atoms, 1.0)
@test from_toml(Hamiltonian, "hamiltonians/h3.toml") == h3

@test to_dict(h1) == OrderedDict{String, Any}(
    "rydberg" => OrderedDict{String, Any}(
        "atoms" => [collect(Int, each) for each in atoms],
        "C" => 1.0,
    ),
    "xterm" => OrderedDict{String, Any}(
        "nsites" => 5,
        "omega" => 0.1,
        "phi" => 0.2,
    ),
    "zterm" => OrderedDict{String, Any}(
        "nsites" => 5,
        "delta" => 0.1,
    )
)

h4 = XTerm([0.1, 0.2, 0.3, 0.4, 0.5], 0.1) + NTerm(5, 0.1) + RydInteract(atoms, 1.0)
@test from_toml(Hamiltonian, "hamiltonians/h4.toml") == h4

atoms = square_lattice(10, 0.8)
space = blockade_subspace(atoms, 1.0)
write_subspace("test.subspace", space)
@test read_subspace("test.subspace") == space
