from matplotlib import projections
import matplotlib.pyplot as plt
import numpy as np

import matplotlib.animation as animation

def plot_trajectory(trajectory, magnets):
    
    ax = plt.figure(figsize=(8, 8)).add_subplot(projection='3d')
    # plt.title(f'Magnetic Pendulum Trajectory')
    ax.set_title(f'Magnetic Pendulum Trajectory')
    # Plot the trajectory
    ax.plot(*trajectory.T, label='Pendulum Path', zorder=1)  

    # Plot the magnet positions
    ax.scatter(*zip(*[m.pos.to_cartesian().coords.astype(np.double) for m in magnets]), color='red', label='Magnets', zorder=2)  

    # Plot the final position of the pendulum
    final_position = trajectory[-1]
    ax.scatter(*final_position, color='cyan', label='Final Position', s=10, zorder=3)

    ax.legend()
    ax.grid(True)
    ax.axis('equal')
    plt.show()

def plot_trajectory_rotate(trajectory, magnets):

    ax = plt.figure(figsize=(8, 8)).add_subplot(projection='3d')
    # plt.title(f'Magnetic Pendulum Trajectory')
    ax.set_title(f'Magnetic Pendulum Trajectory')
    # Plot the trajectory
    ax.plot(*trajectory.T, label='Pendulum Path', zorder=1)  

    # Plot the magnet positions
    ax.scatter(*zip(*[m.pos.to_cartesian().coords.astype(np.double) for m in magnets]), color='red', label='Magnets', zorder=2)  

    # Plot the final position of the pendulum
    final_position = trajectory[-1]
    ax.scatter(*final_position, color='cyan', label='Final Position', s=10, zorder=3)

    ax.legend()
    ax.grid(True)
    ax.axis('equal')


    # Rotate the axes and update
    while True:
        for angle in range(0, 360):
            # Normalize the angle to the range [-180, 180] for display
            angle_norm = (angle + 180) % 360 - 180
            # Cycle through a full rotation of elevation, then azimuth, roll, and all
            elev = azim = roll = 0
            elev = 13
            azim = angle_norm


            # Update the axis view and title
            ax.view_init(elev, azim, roll)
            plt.title('Elevation: %d°, Azimuth: %d°, Roll: %d°' % (elev, azim, roll))

            plt.draw()
            plt.pause(.001)

def plot_trajectory_animation(trajectory, magnets, dt=0.01):

    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(projection="3d")
    # plt.title(f'Magnetic Pendulum Trajectory')
    ax.set_title(f'Magnetic Pendulum Trajectory')
    # Plot the trajectory
    ax.plot(*trajectory.T, label='Pendulum Path', color='gray', alpha=0.1, zorder=1)  

    # Plot the magnet positions
    ax.scatter(*zip(*[m.pos.to_cartesian().coords.astype(np.double) for m in magnets]), color='red', label='Magnets', zorder=2)  

    # Plot the final & initial positions of the pendulum
    final_position = trajectory[-1]
    initial_position = trajectory[0]
    ax.scatter(*final_position, color='cyan', label='Final Position', s=10, zorder=3)
    ax.scatter(*initial_position, color='black', label='Initial Position', s=10, zorder=3)

    ax.legend()
    ax.grid(True)
    ax.axis('equal') 

    line, = ax.plot([], [], [], 'o-', lw=2)
    trace, = ax.plot([], [], [], '.-', lw=1, ms=2)
    time_template = 'time = %.1f'
    time_text = ax.text(0.05, 0.9, 0.9, '', transform=ax.transAxes)


    def animate(i):
        this = trajectory[i:i+1]
        history = trajectory[:i]

        line.set_data_3d(this.T)
        trace.set_data_3d(history.T)
        time_text.set_text(time_template % (i))
        return line, trace, time_text

    N = len(trajectory)
    ani = animation.FuncAnimation(
        fig, animate, N, interval=dt/3, blit=True)
    plt.show()