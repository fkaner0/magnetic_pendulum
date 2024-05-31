import itertools
import numpy as np
from math import sin, cos

def run_simulation(pendulum, dt, max_steps=None):
    trajectory = []  # Initialise blank list

    if max_steps:
        myiter = range(max_steps)
    else:
        myiter = itertools.count()

    for i in myiter:
        pos, theta_dot, phi_dot = pendulum.update(dt)  # Update the pendulum position
        trajectory.append(pos.copy().astype(np.double))  # Store the current position
        speed = np.linalg.norm(get_cartesian_velocity(*pos, theta_dot, phi_dot))
        if speed < 1e-7:
            break  # Stop if pendulum comes to rest

        if (i % 100000/dt) == 0:
            print(f't={i//dt},\tv={speed}')
 
    return np.array(trajectory)  # Return as numpy array

def get_cartesian_velocity(r, theta, phi, theta_dot, phi_dot):
    x_dot = r*(cos(theta)*theta_dot*cos(phi) - sin(theta)*sin(phi)*phi_dot)
    y_dot = r*(cos(theta)*theta_dot*cos(phi) + sin(theta)*cos(phi)*phi_dot)
    z_dot = r*(sin(theta))
    return np.array([x_dot, y_dot, z_dot])