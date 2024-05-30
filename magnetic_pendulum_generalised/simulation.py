import numpy as np

def run_simulation(pendulum, dt, max_steps):
    trajectory = []  # Initialise blank list
    for _ in range(max_steps):
        pos = pendulum.update(dt)  # Update the pendulum position
        trajectory.append(pos.copy().astype(np.double))  # Store the current position
        # if np.linalg.norm(pendulum.vel.to_cartesian().coords[-1]) < 1e-3:
        #     break  # Stop if pendulum comes to rest
    return np.array(trajectory)  # Return as numpy array