# Advection Equation Solver

## Introduction
This code is designed to numerically solve the steady-state
 <a href="https://en.wikipedia.org/wiki/Advection">advection equation</a> 
using the <a href="https://en.wikipedia.org/wiki/Discontinuous_Galerkin_method">
discontinuous Galerkin (DG) method</a>.


**Note**
This is refactored code from the 
<a href="https://www.dealii.org/8.4.0/doxygen/deal.II/step_30.html">
step-30</a> in the deal.ii tutorial. The capabilities of anisotropic 
mesh refinement has been removed to more clearly show how to implement DG methods using
the deal.ii library.


## Requirements
The requirements for this software is 
<a href="https://www.dealii.org">deal.ii library</a> version 8.4.0 or higher,
<a href="https://www.cmake.org">CMake</a> version 2.8 or higher.

## Installation
First obtain and install a copy of the dealii
<a href="https://www.dealii.org">deal.ii library</a> version 8.4.0 or higher. 

## Compiling
To generate a makefile for this code using CMake type into the terminal:

*cmake . -DDEAL_II_DIR=/path_to_deal.ii*

To can compile the code use:

*make release*

## Running
To run the executable use:

*./main*


