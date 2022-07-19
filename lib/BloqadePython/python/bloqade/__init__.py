__version__ = '0.1.0'
import juliacall
import numpy as np
juliacall.Main.seval('using Bloqade')
# juliacall.Main.seval('using BloqadePython')
bloqade = juliacall.Main.seval('Bloqade')

X = bloqade.X
Y = bloqade.Y
Z = bloqade.Z

def rydberg_h(atom_positions, C = 2 * np.pi * 862690, omega=None, phi = None, delta = None):
    return bloqade.rydberg_h(atom_positions, C, omega, phi, delta)

def _convert_pylist_atom_positions(atoms):
    positions = np.array(atoms)
    ndims = len(positions.shape)
    positions = positions.reshape((positions.shape[0], 1))
    return [each for each in positions]

def unit_disk_graph(atoms, radius = 1):
    atoms = _convert_pylist_atom_positions(atoms)
    return bloqade.unit_disk_graph(atoms, radius)

def mis_postprocessing(config, graph, ntrials : int =10):
    return bloqade.mis_postprocessing(config, graph, ntrials=ntrials)



class AbstractBlock(juliacall.AnyValue):

    def __repr__(self) -> str:
        return juliacall.Main.repr(self)
