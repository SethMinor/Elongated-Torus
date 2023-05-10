# Elongated-Torus
Repo for MATLAB material.

## *vortex_simulation.m*
> Driver code for vortex simulation; choose how many vortices you'd like to see the dynamics of here, as well as their starting locations (isothermal coordinates).

## *jacobitheta1.m* and *Djacobitheta1.m*
> Function files to compute the first several terms (double precision is usually reached before a truncation of 15 or so terms) of the first Jacobi theta function and its first derivative. 

## *vortex_velocity.m* and *complexpotential.m*
> MATLAB functions that computes the physical velocity of a set of vortices and their corresponding complex potential function, respectively, using a relation involving Jacobi theta functions is used (see [_Fetter et. al._](https://journals.aps.org/pra/pdf/10.1103/PhysRevA.101.053606?casa_token=Y-7DK7Ny6GYAAAAA%3A6d0WPKGSS2jhegscwXxLSe6u0O6XRoSd-A1o1ET2RzNMRYmkRlpXAkEkiH7Ydck_I-JDhGq016_pfQ) for more details). 

## *hamiltonian.m*
> This function computes the total energy of a system of vortices; it also returns classical, quantum and surface curvature contributions.

## *hamiltonian_contours.m*
> Contours and isosurfaces of vortex dipole Hamiltonian.

## *plotwrapped.m*
> A MATLAB function to plot data wrapped on the surface of a torus nicely (based on MATLAB's `wrapToPi` command).

## *ManyOrbitsOneTorus.m*
> A MATLAB function to plot multiple trajectories simultaneously as a parameter is varied, integrated using *vortex_velocity.m*.

## *fixed_points.m*
> A script that takes an initial guess (seed) of a fixed point `zstar` of the torus ODE's and iteratively searches for a solution to the corresponding nonlinear least squares problem (using the Levenberg-Marquardt algorithm option).
