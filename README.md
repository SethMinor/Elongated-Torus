# Elongated-Torus
Public MATLAB repo for projects about superfluid vortices embedded on an elongated torus.

## *simulator.m*
> Driver code for vortex simulation on the surface of a torus or elongated torus. Choose how many vortices you'd like to see the dynamics of here, as well as their starting locations (isothermal coordinates).

## *jacobitheta1.m* and *Djacobitheta1.m*
> Function files to compute the first several terms (double precision is usually reached before a truncation of 15 or so terms) of the first Jacobi theta function and its first derivative.

## *vortex_velocity.m* and *complexpotential.m*
> MATLAB functions that computes the physical velocity of a set of vortices and their corresponding complex potential function, respectively, using a relation involving Jacobi theta functions is used (see Fetter et al. paper). 

## *hamiltonian.m*
> This function computes the total energy of a system of vortices; it also returns classical, quantum and surface curvature contributions.

## *hamiltonian_contours.m*
> Contours and isosurfaces of vortex dipole Hamiltonian.

## *lyapunov_exponents.m* and *Lyapunov_TEST_Lorenz.m*
> Code for numerically approximating Lyapunov exponents of dipole orbit on a torus. Uses continuous QR factorization with GS orthonormalization (see Geist paper). `Lyapunov_TEST_Lorenz.m` tests this code on the Lorenz '63 system.

## *MLE_continuation.m*
> Numerical continuation code for tracking a branch of orbits as torus parameters are perturbed. 

## *plotwrapped.m*
> A MATLAB function to plot data wrapped on the surface of a torus nicely (based on MATLAB's `wrapToPi` command).

## *fixed_points.m*
> A script that takes an initial guess (seed) of a fixed point `wstar` of the torus ODE's and iteratively searches for a solution to the corresponding nonlinear least squares problem (using the Levenberg-Marquardt algorithm option).
