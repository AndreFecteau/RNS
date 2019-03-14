# RNS

This one-dimensional implicit numerical reactive Navier-Stokes solver was written in the goal of completion of my M.A.Sc in Mechanical Engineering at the university of ottawa. It provides 1-Dimensional steady states solution of deflagration using Finite element centered-difference scheme or Godonuv-type solver using the HLLE Reimann Solver.

The Reactive part of the solver can be removed in the ChaiScript.h file commenting out #SOURCES and probably also #RECENTER_FLAME. It "should" work.

## Getting Started

Clone, have the correct libraries listed in Prerequisites and Go!!

### Prerequisites

Eigen3

Chaiscript V.6.1.0

std=c++17

cereal

Openmp
