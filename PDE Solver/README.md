# PDE Solver

The header file `PDESolver.h` contains three classes: `Advection` `Parabolic` `Hyperbolic`.

To solve one of these three PDEs Boundary condition is needed. The equation can be solved using Dirichlet boundary condition. A class `BoundaryCondition` is created that publicly inherits the class for the equation of interest. The class `BoundaryCondition` creates a function 
```
static float function(float x_vaoue)
```
The equation classes use this function to evaluate the boundary values of the solution.