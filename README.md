# OdeVec

OdeVec is a special-purpose ode-solver. It solves many stiff equation systems
with small- to medium-sized Jacobians optimized on throughput. This is reached through vectorization and
a Python pre-processor that prepares the stiff and sparse equation system symbolically and writes
down the outcome in Fortran syntax in source-file templates.

    Version : 0.1.0
    Author  : Jannes Klee
    Contact : jklee@astrophysik.uni-kiel.de

OdeVec is orginally written to solve chemical networks within the vectorized hydro-code [Fosite](https://github.com/tillense/fosite).

The implementation is based on the algorithm described in
[Description of LSODE](https://computing.llnl.gov/sites/default/files/ODEPACK_pub2_u113855.pdf),
which is included in [ODEPACK](https://computing.llnl.gov/projects/odepack/publications), by A. C. Hindmarsh.
It is a [BDF](https://en.wikipedia.org/wiki/Backward_differentiation_formula)-method.

# General Usage Notes

Run the preprocessor, e.g., with

    ./pre_odevec.py --nvector=32 --example="ROBER" && cd build

This creates the necessary fortran files, with a vector length of 32 and the rober example.

The files need to be compiled by hand at the moment (or included into a larger framework):

    gfortran -O2 -c -g -cpp -fcheck=all -fno-range-check -Wall -fbacktrace odevec.f90 odevec_commons.f90 test.f90

and linking

    gfortran -O2 -o rober *.o

Finally run

    ./primordial

