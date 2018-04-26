# OdeVec

OdeVec is a special-purpose ode-solver. It solves many independent stiff equation systems
with relatively small Jacobians optimized on throughput. This is reached through vectorization and
a Python pre-processor that prepares the stiff and sparse equation system symbolically and writes
down the outcome in Fortran syntax in source-file templates.

    Version : 0.0
    Author  : Jannes Klee
    Contact : jklee@astrophysik.uni-kiel.de

OdeVec is orginally written to solve chemical networks as passives scalars within the hydro-code
Fosite.

Info: Sparsity is not included yet.

# General Usage Notes

Run the preprocessor, e.g., with

    ./pre_odevec.py --nvector=32 --example="PRIMORDIAL"

This creates the necessary fortran files, with a vector length of 32 and a primordial network
example.

This file needs to be compiled by hand at the moment, e.g. for compilation

    gfortran -O2 -c -g -cpp -DHAVE_COEFFICIENTS -fcheck=all -fno-range-check -Wall -fbacktrace odevec.f90 odevec_commons.f90 test.f90

and linking

    gfortran -o primordial *.o

Finally run

    ./primordial

