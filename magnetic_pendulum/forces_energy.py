from typing import Tuple
from coordinate_systems import CurvilinearCoordinateSystem
import sympy as sp
# from numpy import ndarray
import numpy as np

class Force():
    def __init__(self, vector: np.ndarray[sp.Expr], coordsystem: CurvilinearCoordinateSystem) -> None:
        """
        vector: should be 3-tuple in terms of coordsystem's .unit_vectors and .U (scalars)
        """
        ### probably adds checks in for safety - not the focus right now
        self.F = vector
        self.coordsystem = coordsystem


class Energy():
    def __init__(self, scalar: sp.Expr, coordsys: CurvilinearCoordinateSystem) -> None:
        """
        scalar: should be a scalar expression in terms of coordsystem's .U (scalars)
        """
        self.E = scalar
        self.coordsys = coordsys

    def to_force(self):
        # grad(E) w.r.t the basis U
        vector = -np.vectorize(sp.diff)(self.E, self.coordsys.U)/self.coordsys.lame_coefficients
        return Force(vector, self.coordsys)