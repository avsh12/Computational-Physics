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

# PDE Solver

## Introduction
The `PDE Solver` is a program that solves Partial Differential Equations (PDEs) using three different classes: `Advection`, `Parabolic`, and `Hyperbolic`. This README file provides an overview of the program and explains how to use it.

## Classes
The program includes the following classes:

1. `Advection`: This class solves the advection equation.
2. `Parabolic`: This class solves the parabolic equation.
3. `Hyperbolic`: This class solves the hyperbolic equation.

## Boundary Condition
To solve any of the three PDEs, a boundary condition is required. The program uses the Dirichlet boundary condition. To implement this, a class called `BoundaryCondition` is created, which publicly inherits the class for the equation of interest. The `BoundaryCondition` class includes a function called `function(float x_value)`, which is used by the equation classes to evaluate the boundary values of the solution.

## Implementation
The program includes the following implementation details:

1. `advection.cpp`: This file implements the `Advection` class to solve the advection equation with a boundary condition at t = 0.
2. `Matrix M(10, 2, 0.02)`: This line of code creates an object to store the solution as a two-dimensional matrix. The first argument specifies the maximum value of `x`, the second argument specifies the maximum value of `t`, and the third argument specifies the step size of `x`. The step size for `t` is calculated using the value of `s` as Δt = s×Δx.

## Matrix Operations
The program includes the following matrix operations:

1. `UpdateMatrix()`: This function creates a matrix of the required dimensions and sets the values defined in the `Matrix` class.
2. `MatrixOperation MatOp(M.Nx, M.Nt, M.dx, M.dt)`: This line of code creates an object to handle operations on instances of the `Matrix` class. The first two arguments are the ranges of `x` and `t`, starting from zero. The next two arguments are the step sizes of `x` and `t`, respectively.

## Setting Boundary Condition
To specify the boundary condition, the program includes the following steps:

1. Call the function `tBoundary(&BC.u_t0, 0, M.u)`. The first argument takes the function defined in the `BoundaryCondition` class. The second argument takes the `t` index of the matrix (zero means t = 0). The third argument takes the matrix into which the boundary condition is to be written.

## Solving the Equation
To solve the equation, the program includes the following steps:

1. Call the `setMatrix` function of the `MatrixOperation` class. This function uses another function defined in the `Advection` class, `functionDef`, to solve the equation. The `setMatrix` function takes the start and last indices of the variables `x` and `t`, as well as the matrix `M.u`.

## Output
The solution is written to a text file using the `fileOutput` function defined in the `MatrixOperation` class.

