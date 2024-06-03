from typing import Self  # NB: this only works for py > 3.11


class CoordinateSystem:
    def __init__(self, coords) -> None:
        self.coords = coords
        ## not finished
    
    def to(self, coord_system: Self) -> Self:
        """convert to another system - not sure if this is good struct"""
        if coord_system == type(self):
            return self
        else:
            new_system = coord_system.__name__
            function_name = f'_to_{new_system.lower()}'
            try:
                converter = getattr(self, function_name)
            except AttributeError:
                raise NotImplementedError
            return converter()
        

class Cartesian(CoordinateSystem):

    def _to_polar(self):
        pass

class Polar(CoordinateSystem):
    """NB: this polar form has theta defined as angle to NEGATIVE z axis"""

    def _to_cartesian(self):
        pass

############## Thinking: how do I change to include stuff like position vectors

# class Polar(CoordinateSystem):

from sympy import symbols
from sympy.vector import CoordSys3D, express
from sympy import sin, cos, Derivative

C = CoordSys3D('C')  # standard Cartesian system
r, θ, ϕ = symbols('r θ ϕ')

r_ = (r*sin(θ)*cos(φ))*C.i + (r*sin(θ)*sin(φ))*C.j + (r*cos(θ))*C.k
print(Derivative(r_, r))

