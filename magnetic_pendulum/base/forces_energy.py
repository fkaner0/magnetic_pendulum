from typing import Tuple
from .coordinate_systems import CurvilinearCoordinateSystem
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
    
    def __add__(self, other):
        if not isinstance(other, Force):
            raise NotImplementedError("must add to type of Force")
        if self.coordsystem != other.coordsystem:
            raise NotImplementedError("must add to Force with same coordinate system")
        return Force(self.F + other.F, self.coordsystem)


class Energy():
    def __init__(self, scalar: sp.Expr, coordsystem: CurvilinearCoordinateSystem) -> None:
        """
        scalar: should be a scalar expression in terms of coordsystem's .U (scalars)
        """
        self.E = scalar
        self.coordsystem = coordsystem

    def to_force(self):
        # grad(E) w.r.t the basis U
        vector = -np.vectorize(sp.diff)(self.E, self.coordsystem.U)/self.coordsystem.lame_coefficients
        return Force(vector, self.coordsystem)

    def __add__(self, other):
        if not isinstance(other, Energy):
            raise NotImplementedError("must add to type of Energy")
        if self.coordsystem != other.coordsystem:
            raise NotImplementedError("must add to Energy with same coordinate system")
        return Energy(self.E + other.E, self.coordsystem)