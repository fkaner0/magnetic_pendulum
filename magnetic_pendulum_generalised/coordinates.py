from abc import abstractmethod
import numpy as np
from math import sin, cos, pi, sqrt


class Coordinates:

    def __init__(self, coords: np.ndarray) -> None:
        coords = np.array(coords, dtype=np.longdouble)
        if coords.shape != (self._coord_length(),):
            raise TypeError(f"{type(self).__name__} needs {self._coord_length()} components. wrong shape.")
        self.coords = coords

    def get_distance_to(self, other):
        if not isinstance(other, type(self)):
            raise TypeError("Distance must be between coordinates of same type")

        return self._distance_to(other)

    @abstractmethod
    def _distance_to(self, other):
        pass

    @property
    @abstractmethod
    def _coord_length(self) -> int:
        pass



class CartesianCoordinates(Coordinates):

    def __init__(self, coords: np.ndarray) -> None:
        super().__init__(coords)
        self.x, self.y, self.z = coords


    def _coord_length(self):
        return 3

    def _distance_to(self, other):
        return np.linalg.norm(self.coords - other.coords)


class PolarCoordinates(Coordinates):

    def __init__(self, coords: np.ndarray) -> None:
        super().__init__(coords)
        coords = self.coords
        # coords_old = coords.copy()
        coords[1] = coords[1] % (2*pi)
        coords[2] = coords[2] % (2*pi)

        if coords[1] >= pi:
            coords[1] = 2*pi - coords[1]
            coords[2] = coords[2] + pi
        coords[2] = coords[2] % (2*pi)

        # if np.all(coords_old == coords):
        #     print('no change')
        # else:
        #     print('converting:')
        #     for i,(old,new) in enumerate(zip(coords_old, coords)):
        #         if old != new:
        #             print(f'\t[{i}]:\t{old}\t{new}')

        self.coords = coords

        self.r, self.theta, self.phi = coords


    def _coord_length(self):
        return 3

    def _distance_to(self, other):
        r1, t1, p1 = self.coords
        r2, t2, p2 = other.coords

        return sqrt(
            (r1)**2 + (r2)**2
            - 2*(r1)*(r2)*(sin(t1)*sin(t2)*cos(p1-p2) + cos(t1)*cos(t2))
        )
    
    def to_cartesian(self) -> CartesianCoordinates:
        x = self.r * sin(self.theta) * cos(self.phi)
        y = self.r * sin(self.theta) * sin(self.phi)
        z = -self.r * cos(self.theta)
        return CartesianCoordinates([x,y,z])


# c1 = CartesianCoordinates(np.array([1, pi/4, pi/6]))
# c2 = CartesianCoordinates(np.array([3, pi/2, pi/4]))

# print(c1.get_distance_to(c2))