# Polynomial Interpolation
The header file `interpolation.h` contains methods for interpolating a given dataset using polynomial function by using Lagrange and Neville algorithms.

Different algorithms have separate classes for them. There are two classes: `LagrangeInterpolation` and `NevilleInterpolation`.

`data_len` specifies the length of the data to be interpolated. 

The functions `setX` and `setY` take the x and y values of the data.

The interpolated data can be retrived using function `InterpolatingPolynomial` at any given value of x.

`polynomial_interpolation.cpp` is a sample program implementing *Neville* algorithm for interpolating dataset generated using a function.
