from functools import partial
import sympy as sp
import numpy as np
from base import CurvilinearCoordinateSystem, Force, Energy

class PendulumForces:
    def __init__(self, coordsystem: CurvilinearCoordinateSystem, **kwargs) -> None:
        self.coordsystem = coordsystem

        # use aliases for energy and force to get around my bad structure...
        self.Energy = partial(Energy, coordsystem=self.coordsystem)
        self.Force = partial(Force, coordsystem=self.coordsystem)

    def gravitational_potential_energy(self) -> Energy:
        m, g, z0 = sp.symbols('m g z0')
        # where z0 is an arbitrary fixed height (s.t. energy is always +ve)
        z = self.coordsystem.position_vector.dot(self.coordsystem.C.k)
        h = z0 + z
        energy_equation = m*g*h
        return self.Energy(energy_equation, coordsystem=self.coordsystem)

    def magnetic_potential_energy_ASSUMPTION(self, num_magnets) -> Energy:
        namerange = lambda name: f'{name}_mag\,0:{num_magnets}'
        names = ' '.join([namerange(u.name) for u in self.coordsystem.U])

        polarities = np.array(sp.symbols(namerange('\\rho')))
        magnet_positions = np.array(sp.symbols(names)).reshape(-1,num_magnets).T

        magnetic_energy_of = lambda magnet_pos: polarities/(self.coordsystem.get_distance(magnet_pos))**3 * sp.Rational(1,3)
        return np.vectorize(magnetic_energy_of, magnet_positions)
    
    def magnetic_potential_ASSUMPTION_single(self) -> Energy:
        polarity = sp.symbols('\\rho')
        magnet_position = sp.symbols(' '.join([f'{u.name}_mag' for u in self.coordsystem.U]))
        energy_equation = polarity/(self.coordsystem.get_distance(magnet_position))**3 * sp.Rational(1,3)
        return self.Energy(energy_equation)

    def magnetic_potential_ASSUMPTION(self, num_magnets) -> Energy:
        for i in range(num_magnets):
            polarity_i = sp.symbols(f'\\rho_{i}')
            magnet_position_i = sp.symbols(' '.join([f'{u.name}_mag\,{i}' for u in self.coordsystem.U]))
            # self.magnetic_potential_ASSUMPTION_single().subs(zip(magnet_position_i))