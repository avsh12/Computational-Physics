# PDE Solver

The header file `PDESolver.h` contains three classes: `Advection` `Parabolic` `Hyperbolic`.

To solve one of these three PDEs Boundary condition is needed. The equation can be solved using Dirichlet boundary condition. A class `BoundaryCondition` is created that publicly inherits the class for the equation of interest. The class `BoundaryCondition` creates a function 
```
static float function(float x_value)
```
The equation classes use this function to evaluate the boundary values of the solution.

The program `advection.cpp` implements the class to solve the advection with boundary condition at t = 0.
The 
```
Matrix M(10, 2, 0.02)
```
creates an object to store the solution as a two dimensional matrix. The first and second arguments specify the maximum values of `x` and `t`. The third argument specifies the step size of `x`. The step size for `t` is calculated using the value of s as Δt = s×Δx. 

The function `UpdateMatrix()` creates the matrix of the required dimension and sets the values defined in the `Matrix` class.


