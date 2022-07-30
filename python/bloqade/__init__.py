__version__ = '0.1.1'
import juliacall
import numpy as np
juliacall.Main.seval('using Bloqade')
juliacall.Main.seval('using PythonCall')

bloqade = juliacall.Main.seval('Bloqade')
pyconvert = juliacall.Main.seval('pyconvert')

X = bloqade.X
Y = bloqade.Y
Z = bloqade.Z

VectorOfInt = juliacall.Main.seval(r'Vector{Int}')
VectorOfFloat64 = juliacall.Main.seval(r'Vector{Float64}')
VectorOfVectorOfFloat64 = juliacall.Main.seval(r'Vector{Vector{Float64}}')
atom_pos_type = VectorOfVectorOfFloat64

def rydberg_h(atom_positions, C = 2 * np.pi * 862690, omega=None, phi = None, delta = None):
    return bloqade.rydberg_h(atom_positions, C, omega, phi, delta)

def _convert_pylist_atom_positions(atoms):
    positions = np.array(atoms)
    ndims = len(positions.shape)
    if ndims == 1: # list of 1D positions
        positions = positions.reshape((positions.shape[0], 1))

    return pyconvert(atom_pos_type, [[coor for coor in pos] for pos in positions])

def unit_disk_graph(atoms, radius = 1):
    atoms = _convert_pylist_atom_positions(atoms)
    return bloqade.unit_disk_graph(atoms, radius)

def mis_postprocessing(config, graph, ntrials : int =10):
    jl_configs = bloqade.mis_postprocessing(pyconvert(VectorOfInt, config), graph, ntrials=ntrials)
    return np.array(jl_configs)

class AbstractBlock(juliacall.AnyValue):

    def __repr__(self) -> str:
        return juliacall.Main.repr(self)
