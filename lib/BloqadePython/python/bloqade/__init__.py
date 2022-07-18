__version__ = '0.1.0'
import juliacall
import numpy as np
juliacall.Main.seval('using Bloqade')
bloqade = juliacall.Main.seval('Bloqade')

X = bloqade.X
Y = bloqade.Y
Z = bloqade.Z

def rydberg_h(atom_positions, C = 2 * np.pi * 862690, omega=None, phi = None, delta = None):
    return bloqade.rydberg_h(atom_positions, C, omega, phi, delta)

class AbstractBlock(juliacall.AnyValue):

    def __repr__(self) -> str:
        return juliacall.Main.repr(self)
