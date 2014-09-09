C++ code for three bodies integration
=====================================

Integrators
-----------

 * Euler
 * Fourth order Runge-Kutta

Input (initial conditions):
--------------------------

file: init.dat

Distance in AU
Velocity in AU/days
Mass in Solar mass

Output (trajectories):

file: output.dat

How to compile (tested on Mac OS X 10.9):
----------------------------------------

The bash program:

> make three 

The same program with openGL visualization

> make three++

The makefile should work on Linux with little modifications.

How to run:
----------

1. Edit init.dat for initial conditions

2. ./three (three++ for openGL version) 

3. Analyse the output.dat file

