# PDE Solver

The header file `PDESolver.h` contains three classes: `Advection` `Parabolic` `Hyperbolic`.

To solve one of these three PDEs Boundary condition is needed. The equation can be solved using Dirichlet boundary condition. A class `BoundaryCondition` is created that publicly inherits the class for the equation of interest. The class `BoundaryCondition` creates a function 
```
static float function(float x_value)
```
The equation classes use this function to evaluate the boundary values of the solution.

The program `advection.cpp` implements the class to solve the advection equation with boundary condition at t = 0.

`Matrix M(10, 2, 0.02)` creates an object to store the solution as a two dimensional matrix. The first and second arguments specify the maximum values of `x` and `t`. The third argument specifies the step size of `x`. The step size for `t` is calculated using the value of s as Δt = s×Δx. 

The function `UpdateMatrix()` creates the matrix of the required dimension and sets the values defined in the `Matrix` class.

`MatrixOperation MatOp(M.Nx, M.Nt, M.dx, M.dt)` creates an object to handle operations on the instances of `Matrix`. The first two arguments are the ranges of `x` and `t` starting from zero. The next two arguments are the step sizes of `x` and `t` respectively.

Next boundary condition is specified by calling the function `tBoundary(&BC.u_t0, 0, M.u)`. The first argument takes the function defined in the class `BoundaryCondition`. The second argument takes the `t` index of the matrix; zero mean t = 0. The third argument takes the matrix into which boundary condition is to be written.

The equation is solved by calling `setMatrix` function of the class `MatrixOperation`. The function uses another function defined in the `Advection` class, `functionDef`, to solve the equation. `setMatrix` takes start and last indices of variables `x` and `t`, and the matrix `M.u`

The solution is written to a text file using the function `fileOutput` defined in `MatrixOperation`.

