
from math import pi, atan
from magnetic_pendulum import MagneticPendulum, run_simulation, plot_trajectory
import numpy as np

from scipy.constants import g

from magnetic_pendulum_generalised import (MagneticPendulum as MagneticPendulum_Polar,
                                           run_simulation as run_simulation_Polar,
                                           plot as plot_Polar)
from magnetic_pendulum_generalised.coordinates import PolarCoordinates
from magnetic_pendulum_generalised.pendulum import Magnet

# Define the magnets in a square centered at the origin

b = 0.3  # Damping coefficient
h = 0.5  # Height above the x-y plane

L = g  # bob length
R = 0.5  # magnet radius
magnets = [(0, R), (0, -R), (-R, 0), (R, 0)]

mag_theta =atan(R/(L+h))
magnets_polar = [Magnet(PolarCoordinates([(L+h), mag_theta, (pi/2) * i]), strength=1) for i in range(4)]


initial_pos_polar = PolarCoordinates([L, 0.2, 0.4])
initial_pos = initial_pos_polar.to_cartesian().coords[:-1]
initial_vel_polar = PolarCoordinates([0,0,0])

# initial_pos = [0.9, -0.6]  # Random fixed initial position
initial_vel = (0, 0)  # Initial velocity
print(f"Initial position: {initial_pos_polar.to_cartesian()}")

# Create the pendulum
pendulum = MagneticPendulum(magnets, b, h, initial_pos, initial_vel)
pendulum_polar = MagneticPendulum_Polar(magnets_polar, b, h, initial_pos_polar, initial_vel_polar)

dt = 0.01  # Time step
max_steps = 10000  # Maximum number of steps
# trajectory = run_simulation(pendulum, dt, max_steps)
trajectory_polar = run_simulation_Polar(pendulum_polar, dt, max_steps)

# plot_trajectory(trajectory, magnets)
# cont = input("any key to continue")

# plot_Polar.plot_trajectory(trajectory_polar, magnets_polar)
# plot_Polar.plot_trajectory_rotate(trajectory_polar, magnets_polar)
# plot_Polar.plot_trajectory_animation(trajectory_polar, magnets_polar, dt)


plot_Polar.plot_trajectory(trajectory_polar, magnets_polar)
plot_Polar.plot_trajectory_rotate(trajectory_polar, magnets_polar)
print("dunzo")