# Elongated-Torus
Public MATLAB repo for projects about superfluid vortices embedded on an elongated torus. Includes both **driver code** and **helper functions**.

Note: download the entire contents of the `Elongated-Torus.zip` file need to ensure a given script works properly.

# Driver code
The fun scripts to run. These do real things.

## **`simulator.m`**
> Driver code for vortex simulation on the surface of a torus or elongated torus. Choose how many vortices you'd like to see the dynamics of here, as well as their starting locations (isothermal coordinates).

## **`PDE_vs_ODE.m`**
> Compares a chosen ODE orbit with the corresponding integrated NLS solution.

## **`hamiltonian_v2.m`** and **`hamiltonian_contours.m`**
> This function computes the total energy of a system of vortices; it also returns classical, quantum and surface curvature contributions. `hamiltonian_contours.m` plots the contours and isosurfaces of vortex dipole Hamiltonians (helpful for selecting constant energy orbits for Poincare sections).

## **`lyapunov_exponents.m`** and **`Lyapunov_TEST_Lorenz.m`**
> Code for numerically approximating Lyapunov exponents of dipole orbit on a torus. Uses continuous QR factorization with GS orthonormalization (see Geist paper). `Lyapunov_TEST_Lorenz.m` tests this code on the Lorenz '63 system.

## **`fixed_points.m`**
> A script that takes an initial fixed point seed `W` and iteratively searches for a fixed point solution `Wstar` to the corresponding nonlinear least squares problem (using the Levenberg-Marquardt algorithm). Also returns eigenvalues and eigenvectors.

## **`MLE_continuation.m`**
> Numerical continuation code for tracking a branch of orbits as torus parameters are perturbed. Also plots the MLE vs eccentricity (epsilon) at different aspect ratios (alpha).

## **`poincare_sections.m`**
> Code for computing Poincare sections on the elongated torus, with an ensemble of initial conditions.

# Helper functions
Likely no need to edit these files.

## **`vortex_velocity_v2.m`** and **`vortex_velocity_lyapunov.m`**
> Function that computes the physical velocity of a set of vortices. `vortex_velocity_lyapunov.m` has special additions for use in computing Lyapunov exponents.

## **`myjacobian.m`**
> Numerically computed Jacobian matrix for ODE. Can use either a first, second, or fourth order derivative approximation. Step length is `epsilon`.

## **`jacobitheta1.m`** and **`Djacobitheta1.m`**
> Helper function files to compute the first several terms (double precision is usually reached before a truncation of 15 or so terms) of the first Jacobi theta function and its first derivative.

## **`plotwrapped.m`**
> Helper function to plot data wrapped on the surface of a torus nicely (based on MATLAB's `wrapToPi` command).

## **`UVwrapped.m`**
> Helper function for wrapping periodic data.
