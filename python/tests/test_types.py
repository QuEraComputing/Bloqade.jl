import bloqade.types as types
from bloqade.types import JLArray, JLMatrix, JLInt32
import numpy as np
import juliacall

def test_JLMatrix_convert():
    assert type(JLMatrix[JLInt32].convert(np.array([[1, 2], [2, 1]]))) == juliacall.ArrayValue
