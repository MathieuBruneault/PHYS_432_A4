# PHYS 432_: Assignment 4

* Name: Mathieu Bruneault
* Python version: 3.x
* The two problems I submit are the advection equation and the hydro solver (whose scripts are respectively name advection_equation.py and hydro_solver.py)
* Collaborators: Alexandre Khoury, Ronan Legin

## List of Files and what they do

1. advection_equation.py

This script numerically solves the advection equation with fixed boundary conditions and an initial distribution f(x) = x. It attempts to solve using two different methods: FTCS and Lax-Friedrich, and displays the progressions in time of these two solutions side by side. From terminal you may run it simply as python advection_equation.py.

2. hydro_solver.py

This script numerically solves the continuity equation and the Euler equation for the 1D motion of sound waves in a uniform density gas. The script plots the density and the velocity field in time. It can be run usong python hydro_solver.py

## Observations/Answers to question

1. Advection Equation

As should be expected, the FTCS scheme is not stable at all, it blow up quite quickly. The Lax-Friedrich method is more well behaved, and it tends to some stationary function as t progresses

2. Hydro solver

With a large enough amplitude and reflective boundary conditions, two waves eventually form at the extremities of our x space, and they then collide at the center, forming a shock. The amplitude of this shock is set by the initial amplitude of density, but the width is most likely set by the initial velocity distribution.
