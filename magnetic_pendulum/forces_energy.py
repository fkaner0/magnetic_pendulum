from typing import Tuple
from coordinate_systems import CurvilinearCoordinateSystem
import sympy as sp\

class Force():
    def __init__(self, vector: Tuple[sp.Expr], coordsystem: CurvilinearCoordinateSystem) -> None:
        """
        vector should be 3-tuple in terms of coordsystem's .unit_vectors and .U (scalars)
        """
        self.F = vector
        self.coordsystem = coordsystem


class Energy():
    def __init__(self, scalar: sp.Expr, coordsys: CurvilinearCoordinateSystem) -> None:
        """
        should be scalar value in terms of coordsystem's .U (scalars)
        """
        self.E = scalar
        self.coordsys = coordsys

    def to_force(self):

        ### haven't really thought through this yet
        vector = []
        for u in self.coordsys.U:
            vector.append(sp.diff(self.E, u))

        return Force(vector, self.coordsys)