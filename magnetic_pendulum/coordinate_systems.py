import sympy as sp
from sympy import sin, cos
from sympy.vector import CoordSys3D, VectorAdd

from functools import cached_property, partial
from typing import Self  # NB: this only works for py > 3.11


class CurvilinearCoordinateSystem:
    def __init__(self, t, symbols, cartesian_position_formula, C: CoordSys3D=CoordSys3D('C'), name: str="") -> None:
        """
        t: sp.Symbol which defines the dynamics of our system
        symbols: 3-tuple of sp.Symbol which we use to define position (these will be converted into functions of t)
        position_vector: 3-tuple, giving the cartesian coordinates of a in the system, in terms of `symbols`
        (optional) C: standard cartesian coordinate system we define this in terms of
        """
        self.name = name
        self.C = CoordSys3D('C')  # standard Cartesian system
        self.t = t
        self.old_symbols = symbols

        # make sure our symbols become (undefined) functions of t
        def func_of_t(*args, **kwargs):
            return sp.Function(*args, **kwargs)(t)
        # tuple of symbols as functions of t
        self.U = sp.symbols(' '.join(s.name for s in symbols), cls=func_of_t)

        position_vector: VectorAdd = (
            cartesian_position_formula[0]*C.i
            + cartesian_position_formula[1]*C.j
            + cartesian_position_formula[2]*C.k
        )
        self.position_vector = position_vector.subs(zip(self.U, self.old_symbols))

    @cached_property
    def unit_vectors(self):
        return [sp.simplify(sp.diff(self.position_vector, u).normalize()) for u in self.U]
        
    @cached_property
    def acceleration_vector(self):
        """
        list giving acceleration in directions of each unit vector
        (note: this distinction becomes clear in spherical coordinates, for example)
        """
        acceleration_formula = sp.diff(self.position_vector, self.t, 2)
        a_in_direction = lambda i: sp.simplify(acceleration_formula.dot(self.unit_vectors[i]))
        return [a_in_direction(i) for i in range(3)]
    
    def acceleration_components(self, force):
        """
        gives acceleration of each of `self.U`, given some force
        force: 3-tuple giving the scalar (???) force in the direction of each unit_vector
        """
        # F = {u: sp.symbols(f'F_{u.name}') for u in self.U}
        U_acceleration = {}
        for i, u in enumerate(self.U):
            # get acceleration in terms of u
            u_a_solns = sp.solve(sp.Eq(self.acceleration_vector[i], force[i]), sp.diff(u, self.t, 2))
            if len(u_a_solns) != 1:
                raise ValueError(f'could not find acceleration in direction of {u} in terms of {sp.diff(u, self.t, 2)}.\nreturned solns: {u_a_solns}')
            U_acceleration[u] = u_a_solns[0]

        return U_acceleration


def _get_cartesian_system():
    symbols = (x, y, z) = sp.symbols('r θ φ')
    position_vector = (x, y, z)

    return partial(CurvilinearCoordinateSystem(
        symbols=symbols,
        cartesian_position_formula=position_vector,
        name='Cartesian'
    ))

def _get_cylindrical_system():
    symbols = (r, θ, z) = sp.symbols('r θ φ')
    position_vector = ((r*cos(θ)), (r*sin(θ)), (z))

    return partial(CurvilinearCoordinateSystem(
        symbols=symbols,
        cartesian_position_formula=position_vector,
        name='Cylindrical'
    ))

def _get_spherical_system():
    symbols = (r, θ, φ) = sp.symbols('r θ φ')
    position_vector = ((r*sin(θ)*cos(φ)), (r*sin(θ)*sin(φ)), (r*cos(θ)))

    return partial(CurvilinearCoordinateSystem(
        symbols=symbols,
        cartesian_position_formula=position_vector,
        name='Spherical'
    ))

def _get_inverted_spherical_system():
    """standard spherical, but θ corresponds to *negative* z-axis"""
    symbols = (r, θ, φ) = sp.symbols('r θ φ')
    position_vector = ((r*sin(θ)*cos(φ)), (r*sin(θ)*sin(φ)), (-r*cos(θ)))

    return partial(CurvilinearCoordinateSystem(
        symbols=symbols,
        cartesian_position_formula=position_vector,
        name='InvSpherical'
    ))


def get_predefined_system(name):
    """[not well implemented, but does the job]
    get a coordinate system from predefined set
    """
    if name == 'Cartesian':
        return _get_cartesian_system()
    elif name == 'Cylindrical':
        return _get_cylindrical_system()
    elif name == 'Spherical':
        return _get_spherical_system()
    elif name == 'InvSpherical':
        return _get_inverted_spherical_system()
    else:
        raise ValueError(f"no predefined system with name '{name}'")
    