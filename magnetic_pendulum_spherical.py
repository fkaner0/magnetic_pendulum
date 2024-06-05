import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp


# PARAMETERS
L = 1.0  # Length of pendulum
b = 0.05  # Damping coefficient
g = 9.81  # Acceleration due to gravity
k = 5.0  # Magnetic strength coefficient
h = 0.5  # Height of the pendulum above the xy-plane
d = 1.0  # Radius of magnets from the origin on the xy-plane
m = 1.0  # Mass of pendulum bob


# Initial conditions
theta0 = np.pi / 4  # 45 degrees
phi0 = np.pi / 4  # 45 degrees
theta_dot0 = 0.0  # Initially stationary
phi_dot0 = 0.0  # Initially stationary
initial_conditions = [theta0, phi0, theta_dot0, phi_dot0]


dt = 0.01  # Time step
t_max = 20  # Maximum time
t_vals = np.arange(0, t_max, dt)


magnet_positions = [
    (d, 0, -(L + h)),
    (-d, 0, -(L + h)),
    (0, -d, -(L + h)),
    (0, d, -(L + h))
]


def magnet_dist(theta, phi, magnet_pos):
    x_m, y_m, z_m = magnet_pos
    x_p = L * np.sin(theta) * np.cos(phi)
    y_p = L * np.sin(theta) * np.sin(phi)
    z_p = -L * np.cos(theta)
    return np.sqrt((x_p - x_m)**2 + (y_p - y_m)**2 + (z_p - z_m)**2)


def mag_potential_energy(theta, phi):
    V = 0
    for pos in magnet_positions:
        r = magnet_dist(theta, phi, pos)
        V += k / r**4  # V_mag is proprtional to 1/dist^4
    return -V


def dVmag_dtheta(theta, phi):
    # Partial derivative of Vmag with respect to theta
    dV_dtheta = 0
    for magnet_pos in magnet_positions:
        r_i = magnet_dist(theta, phi, magnet_pos)
        x_m, y_m, z_m = magnet_pos
        x_p = L * np.sin(theta) * np.cos(phi)
        y_p = L * np.sin(theta) * np.sin(phi)
        z_p = -L * np.cos(theta)
        dri_d_theta = (L * np.cos(theta) * np.cos(phi) * (x_p - x_m) +
                         L * np.cos(theta) * np.sin(phi) * (y_p - y_m) +
                         L * np.sin(theta) * (-z_p - z_m)) / r_i
        dV_dtheta -= 4 * k * dri_d_theta / r_i**5
    return dV_dtheta


def dVmag_dphi(theta, phi):
    # Partial derivative of Vmag with respect to phi
    dV_dphi = 0
    for magnet_pos in magnet_positions:
        r_i = magnet_dist(theta, phi, magnet_pos)
        x_m, y_m, z_m = magnet_pos
        x_p = L * np.sin(theta) * np.cos(phi)
        y_p = L * np.sin(theta) * np.sin(phi)
        z_p = -L * np.cos(theta)
        d_r_i_d_phi = (-L * np.sin(theta) * np.sin(phi) * (x_p - x_m) +
                       L * np.sin(theta) * np.cos(phi) * (y_p - y_m)) / r_i
        dV_dphi -= 4 * k * d_r_i_d_phi / r_i**5
    return dV_dphi


def derivatives(t, y):
    theta, phi, theta_dot, phi_dot = y

    # Partial derivatives of the magnetic potential energy
    dV_dtheta = dVmag_dtheta(theta, phi)
    dV_dphi = dVmag_dphi(theta, phi)

    # Euler Lagrange equations, isolating second derivatives
    theta_double_dot = (m * L**2 * np.sin(theta) * np.cos(theta) * phi_dot**2
                        - m * g * L * np.sin(theta)
                        - dV_dtheta
                        - b * theta_dot) / (m * L**2)
    phi_double_dot = (-2 * m * L**2 * np.sin(theta) * np.cos(theta) * theta_dot * phi_dot
                      - dV_dphi
                      - b * phi_dot) / (m * L**2 * np.sin(theta)**2)
    return [theta_dot, phi_dot, theta_double_dot, phi_double_dot]


# Numerical solution
solution = solve_ivp(derivatives, [0, t_max], initial_conditions, t_vals=t_vals, method='RK45')
theta = solution.y[0]
phi = solution.y[1]


# Convert spherical coordinates to Cartesian coordinates for plotting
x = L * np.sin(theta) * np.cos(phi)
y = L * np.sin(theta) * np.sin(phi)
z = -L * np.cos(theta)


# Plotting the trajectory
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(x, y, z, label='Trajectory')
ax.scatter(x[-1], y[-1], z[-1], color='orange', s=50, label='Final Position', zorder=5)


# Magnet positions
for magnet_pos in magnet_positions:
    ax.scatter(*magnet_pos, color='red', s=50, label='Magnet' if 'Magnet' not in ax.get_legend_handles_labels()[1] else "", zorder=4)


# Labels and legend
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.legend()
plt.show()