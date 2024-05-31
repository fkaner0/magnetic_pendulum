# from abc import abstractmethod
from typing import List
import numpy as np

from magnetic_pendulum_generalised.coordinates import Coordinates, PolarCoordinates

from math import pi, sin,cos
from scipy.constants import g



class Magnet:
    def __init__(self, position:Coordinates, polarity=1, strength=1):
        self.pos = position
        self.polarity = polarity
        self.strength = strength
        


class MagneticPendulum:

    coord_type = PolarCoordinates

    def __init__(self, magnets:List[Magnet], b, h, initial_pos, initial_vel):
        # self.magnets = np.array(magnets, dtype=float)  # Magnet positions in numpy array
        self.magnets = magnets
        self.b = b  # Damping constant
        self.h = h  # Distance of bob above x-y plane

        if not isinstance(initial_pos, self.coord_type):
            initial_pos = self.coord_type(initial_pos)
        self.pos = initial_pos  # Initial position (x, y)

        if not isinstance(initial_vel, self.coord_type):
            initial_vel = self.coord_type(initial_vel)
        # self.vel = initial_vel  # Initial position (x, y)
        self.thetadot = initial_vel.theta
        self.phidot = initial_vel.phi
        self.old_acceleration = np.array([0,0,0], dtype=np.longdouble)

        self.m = 1
        # self.strength = 5
        # self.bob = Magnet(strength, polarity=1)


    # def _distance(self, magnet_pos):
    #     """Magnet distances : first derived equation."""
    #     return self.pos.get_distance_to(magnet_pos)
        # return np.sqrt((magnet_pos[0] - self.pos[0])**2 + (magnet_pos[1] - self.pos[1])**2 + self.h**2)

    def _magnetic_force(self, magnet:Magnet):
        """Magnetic force vector : second derived equation."""
        dist = self.pos.get_distance_to(magnet.pos)
        multiplier = magnet.strength * magnet.polarity * (1/dist**5)

        th1 = self.pos.theta
        th2 = magnet.pos.theta
        ph1 = self.pos.phi
        ph2 = magnet.pos.phi
        dtheta = magnet.pos.r*(cos(th1)*sin(th2)*cos(ph1-ph2) - sin(th1)*cos(th2))
        dphi   = -magnet.pos.r*(sin(th2)*sin(ph1-ph2))

        force = multiplier * np.array([0, dtheta, dphi])

        return force

        # return self.strength * (magnet_pos - self.pos) * (1 / dist**5)

    def _gravitational_force(self):
        """Gravitational force on bob."""
        return self.m * g * np.array([0, -sin(self.pos.theta), 0])

    def _damping_force(self):
        """Damping force on bob."""
        return -self.b * self.pos.r * np.array([0, self.thetadot, self.phidot * self.pos.theta])

    def _total_force(self):
        """Total force acting on bob."""
        magnetic_forces = np.sum([self._magnetic_force(m) for m in self.magnets], axis=0)
        # return self._gravitational_force() + magnetic_forces
        return self._gravitational_force() + self._damping_force() + magnetic_forces


    def _acceleration_components(self, dt):
        

        force = self._total_force()
        L = self.pos.r
        theta = self.pos.theta
        theta_dot = self.thetadot
        phi_dot = self.phidot

        # for theta tiny we don't consider acceleration - everything becomes undefined
        if theta < (0.00000005 * dt):
            return self.old_acceleration

        theta_dotdot = (force[1]/self.m + L*sin(theta)*cos(theta)*((phi_dot)**2))/L
        phi_dotdot = (force[2]/self.m - 2*L*theta_dot*phi_dot*cos(theta))/(L*sin(theta))
        
        # print(f'F  =\t{force}')
        # print(f'ph..:\t{phi_dotdot}\t{np.isfinite(phi_dotdot)}')
        # print(f'ph.\t{phi_dot}')
        # print(f'th..:\t{theta_dotdot}\t{np.isfinite(theta_dotdot)}')
        # print(f'th.\t{theta_dot}')
        # print('\n\n')
        if not np.isfinite(phi_dotdot):
            print('aa')
            # raise ValueError
        

        acceleration = np.array([0, theta_dotdot, phi_dotdot])
        return acceleration

    def update(self, dt):
        """At each time interval, record position and velocity of bob using Euler's Method."""
        # if self.pos.theta < 0.3:
        # print(f'{self.thetadot}\t{self.phidot}')
        # print(f'\t{self.pos.coords}\n\n')
        # if not (np.isfinite(self.thetadot) and np.isfinite(self.phidot)):
        #     print("aaaaaaaaaaa")
        #     raise ValueError
        
        acceleration = self._acceleration_components(dt)
        self.old_acceleration = acceleration
        # print(acceleration)
        vel = np.array([0, self.thetadot, self.phidot])
        vel += acceleration * dt
        _, self.thetadot, self.phidot = vel
        self.pos = self.coord_type(self.pos.coords + (vel * dt))
        return (self.pos.to_cartesian().coords, self.thetadot, self.phidot) 