# Zerilli

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://svretina.github.io/Zerilli.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://svretina.github.io/Zerilli.jl/dev/)
[![Build Status](https://github.com/svretina/Zerilli.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/svretina/Zerilli.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![Coverage](https://codecov.io/gh/svretina/Zerilli.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/svretina/Zerilli.jl)

This package is evolving the radial part of the Zerilli equation to calculate 
the evolution of black hole pertubations. The equations are put into first order
formulation and boundary conditions are applied via the characteristic approach. Boundary conditions
can be either reflective or radiative. It uses 2nd order stencils for the discretization
in the interior and 1st order at the boundary points, so that we employ stability proofs from 
the Summation By Parts (SBP) property. 


## TODO
- write some docs 
- write correct initial data
- write potential
- put correct equations with the potential term
- (add threading? not needed right now)
