"""
1 dimensional hydro solver using the donor cell advection scheme

@author: Mathieu Bruneault
Feb. 24th 2020
"""

import numpy as np
import matplotlib.pyplot as plt

# Initiate some parameters
del_x = 0.001
del_t = 0.000001
n = 10000
time_it = 100000
cs = 300

def gaussian(x, mu, sig):
    # Returns a Gaussian distribution with mean mu and standard deviation sig, evaluated at x
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

# Initiate f1 and f2 as Gaussian distributions
x = np.linspace(0, n * del_x, num = n)
f1 = 0.5*gaussian(x, np.mean(x), n * del_x / 10) + 0.5
f2 = 100 * (f1 - 0.5)

# Set up the animation and draw the initial fuctions
plt.ion()
fig, axes = plt.subplots(1,2)
axes[0].set_title('Density')
axes[1].set_title('Velocity Field')
axes[0].set_ylim([0,3])
axes[1].set_ylim([-300,300])
x1, = axes[0].plot(x, f1)
x2, = axes[1].plot(x, f2 / f1)

def get_current(u_bound, f):
    # Find the boundaries at which u<0 vs u>0 and compute the Js appropriately depending on the sign of u
    J = np.zeros(n-1)
    neg_u = np.array(np.where(u_bound < 0))
    pos_u = np.array(np.where(u_bound > 0))
    J[pos_u] = u_bound[pos_u] * f[pos_u]
    J[neg_u] = u_bound[neg_u] * f[neg_u + 1]
    return J

def update_functions(f1, f2, J1, J2):
    # Update each cell not at the boundary based on the currents at its two extremities
    f1[1:n-1] -= del_t / del_x * (J1[1 : n-1] - J1[0 : n - 2])
    f2[1:n-1] -= del_t / del_x * (J2[1 : n-1] - J2[0 : n - 2])
    f2[1:n-1] -= del_t / del_x * cs ** 2 * (f1[2:] - f1[:n-2])

    # Apply reflective boundary conditions to both f1 and f2
    f1[0] -= (del_t / del_x) * J1[0]
    f1[-1] += (del_t / del_x) * J1[-1]
    f2[0] -= (del_t / del_x) * J2[0]
    f2[-1] += (del_t / del_x) * J2[-1]

def draw(f1, f2):
    # Update the plot with new functions
    x1.set_ydata(f1)
    x2.set_ydata(f2 / f1)
    fig.canvas.draw()
    plt.pause(0.001)

for i in np.arange(time_it):
    # Compute the velocity at all boundaries between two cells (not at the extremities)
    u_bound = 0.5 * (f2[0 : n-1] / f1[0 : n-1] + f2[1 : n] / f1[1 : n])

    # Get the currents and update the functions
    J1 = get_current(u_bound,f1)
    J2 = get_current(u_bound,f2)
    update_functions(f1, f2, J1, J2)

    # Draw the plots every 100 iterations (I found this looks nicest)
    if i%100 == 0:
        draw(f1, f2)
