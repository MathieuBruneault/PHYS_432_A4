"""
Forward-Time Central Space and Lax Friedrichs methods to solve the advection equation

@author: Mathieu Bruneault
Feb. 24th 2020
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Initiate some parameters
del_x = 0.1
u = -0.1
num_cells = 100
time_it = 2000

# Initiate the function f we are solving for at the first time iteration, by setting it equal to x
x = np.linspace(0, num_cells * del_x, num = num_cells)
f0 = x.copy()

# These lists keep track of the system at each iteration, and are an artefact of me saving the
# animation rather than displaying it
f_list_ftcs = np.array([f0])
f_list_lax = np.array([f0])

# Compute the time step and a parameter alpha that will make updates easier on the eyes later in the code
# I use a smaller time step than the biggest permitted for stability, meaning it will also be stable
del_t = del_x / (10 * np.abs(u))
alpha = (u * del_t / (2*del_x))

# Set up the animation and draw the initial fuctions
plt.ion()
fig, axes = plt.subplots(1,2)
axes[0].set_title('FTCS')
axes[1].set_title('Lax-Friedrich')
x1, = axes[0].plot(x, f0)
x2, = axes[1].plot(x, f0)
fig.canvas.draw()


def ftcs_update(x, f):
    # Update the function through one time step using the FTCS scheme. This uses np.roll, which will
    # cause periodic boundary conditions, so the bondaries need to be updated
    f -= alpha * (np.roll(f, -1) - np.roll(f, 1))
    f[0] = x[0]
    f[-1] = x[-1]
    return f

def lax_update(x, f):
    # Update the function through one time step using the Lax-Friedrich scheme. This uses np.roll, which will
    # cause periodic boundary conditions, so the bondaries need to be updated
    f = 1 / 2 * (np.roll(f, -1) + np.roll(f, 1)) - alpha * (np.roll(f, -1) - np.roll(f, 1))
    f[0] = x[0]
    f[-1] = x[-1]
    return f

def draw(f1, f2):
    # Update the plot with new functions
    x1.set_ydata(f1)
    x2.set_ydata(f2)
    fig.canvas.draw()
    plt.pause(0.001)

# Loop through time iterations
for i in np.arange(time_it):
    # Update f1 and f2 and append them to their respective lists
    f1 = ftcs_update(x, f_list_ftcs[-1])
    f_list_ftcs = np.append(f_list_ftcs, [f1], axis = 0)
    f2 = lax_update(x, f_list_lax[-1])
    f_list_lax = np.append(f_list_lax, [f2], axis = 0)

    # Draw the plots every 5 iterations (I found this looks nicest)
    if i % 5 == 1:
        draw(f1, f2)
