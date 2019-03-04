#!/usr/bin/env python
""" Python preprocessor for the stiff BDF solver OdeVec

This script sets up the system to solve and does some precomputations
if possible. This includes writing down Fortran code in *.F90 template
files. There, placeholder (Pragmas) are placed, which are searched and
replaced. The precomputations includes the symbolic computation of the
Jacobian and preordering of the latter.
"""

from __future__ import print_function
from sympy import symbols, Matrix, eye, zeros, fcode, SparseMatrix, lambdify, nsimplify
from scipy import sparse
from shutil import copyfile
from subprocess import check_output
import os
import numpy as np
import argparse


def GetSystemToSolve(nvector,example,kromefile=""):
    """ Sets up the right-hand-side of the problem

    Defines the RHS that should be solved. Currently there are three
    possibilties:
    1. Robertson's example - 3 chemical species (very stiff)
    2. A primordial network - taken from the KROME package
    3. Arbitrary networks with KROME

    Parameters
    ----------
    nvector : int
        The used vector length.
    example : str
	The chosen examples "ROBER", "PRIMORDIAL" or "KROME"
    kromefile : str, optional
	Path to kromefile were the RHS of the network is dumped in
	Python syntax (default is empty string)

    Returns
    -------
    y : sympy matrix
        The functions of the problem.
    rhs : sympy matrix
        The right-hand side of the problem.
    nvector : int
        Vector length that is intended to use.
    neq : int
        Number of equations. Length of array y and rhs.
    maxorder : int
        Maximum order until LU decomposition should be done
        in precompution.
    nrea : int
        Number of reactions. Only necessary for primordial setups.
    """

    # general global constants and printouts
    nvector = nvector
    maxorder = 5
    print("#  Setting up system...")
    print("#  Vector length \"nvector\":", nvector)

    #------------------------------------------------------------------------------#
    # Robertson's example
    if (example=="ROBER"):
        print("#  Test: ROBER")

        neq = 3
        k = 0
        nrea = 0

	# define sympy symbols
        y = symbols('y(i\,1:%d)'%(neq+1))

        rhs = Matrix([-0.04*y[0]+1e4*y[1]*y[2],
		      0.04*y[0]-3e7*y[1]*y[1]-1e4*y[1]*y[2],
                      3e7*y[1]*y[1]])

    # Primordial network with 16 different species
    elif (example=="PRIMORDIAL"):
        print("#  Test: PRIMORDIAL")

        neq = 16
        nrea = 38

	# indexing of species
        idx_E = 0
        idx_Hk = 1
        idx_H = 2
        idx_HE = 3
        idx_H2 = 4
        idx_D = 5
        idx_HD = 6
        idx_Hj = 7
        idx_HEj = 8
        idx_H2j = 9
        idx_Dj = 10
        idx_HEjj = 11
        idx_CR = 12
        idx_g = 13
        idx_Tgas = 14
        idx_dummy = 15

	# define sympy symbols
        y = symbols('y(i\,1:%d)' % (neq + 1))
        k = symbols('k(i\,1:%d)' % (nrea + 1))

        # RHS of primordial network
        rhs = Matrix([- k[0] * y[idx_H] * y[idx_E] + 2.e0 * k[0] * y[idx_H] * y[idx_E]
                      - k[1] * y[idx_Hj] * y[idx_E] - k[2] * y[idx_Hj] * y[idx_E]
                      - k[3] * y[idx_HE] * y[idx_E] + 2.e0 * k[3] * y[idx_HE] * y[idx_E]
                      - k[4] * y[idx_HEj] * y[idx_E] - k[5] * y[idx_HEj] * y[idx_E]
                      - k[6] * y[idx_HEj] * y[idx_E] + 2.e0 * k[6] * y[idx_HEj] * y[idx_E]
                      - k[7] * y[idx_HEjj] * y[idx_E]
                      - k[8] * y[idx_H] * y[idx_E]
                      + k[9] * y[idx_Hk] * y[idx_H]
                      + k[10] * y[idx_Hk] * y[idx_H]
                      - k[15] * y[idx_H2] * y[idx_E] + k[15] * y[idx_H2] * y[idx_E]
                      - k[17] * y[idx_Hk] * y[idx_E] + 2.e0 * k[17] * y[idx_Hk] * y[idx_E]
                      + k[18] * y[idx_Hk] * y[idx_H]
                      + k[19] * y[idx_Hk] * y[idx_H]
                      + k[21] * y[idx_Hk] * y[idx_Hj]
                      - k[22] * y[idx_H2j] * y[idx_E]
                      - k[23] * y[idx_H2j] * y[idx_E]
                      + k[36] * y[idx_D] * y[idx_Hk]
                      - k[37] * y[idx_Dj] * y[idx_E],
                      + k[8] * y[idx_H] * y[idx_E]
                      - k[9] * y[idx_Hk] * y[idx_H]
                      - k[10] * y[idx_Hk] * y[idx_H]
                      - k[17] * y[idx_Hk] * y[idx_E]
                      - k[18] * y[idx_Hk] * y[idx_H]
                      - k[19] * y[idx_Hk] * y[idx_H]
                      - k[20] * y[idx_Hk] * y[idx_Hj]
                      - k[21] * y[idx_Hk] * y[idx_Hj]
                      - k[24] * y[idx_H2j] * y[idx_Hk]
                      - k[36] * y[idx_D] * y[idx_Hk],
                      - k[0] * y[idx_H] * y[idx_E]
                      + k[1] * y[idx_Hj] * y[idx_E]
                      + k[2] * y[idx_Hj] * y[idx_E]
                      - k[8] * y[idx_H] * y[idx_E]
                      - k[9] * y[idx_Hk] * y[idx_H]
                      - k[10] * y[idx_Hk] * y[idx_H]
                      - k[11] * y[idx_H] * y[idx_Hj]
                      - k[12] * y[idx_H] * y[idx_Hj]
                      - k[13] * y[idx_H2j] * y[idx_H]
                      + k[14] * y[idx_H2] * y[idx_Hj]
                      + 2.e0 * k[15] * y[idx_H2] * y[idx_E]
                      - k[16] * y[idx_H2] * y[idx_H] + 3.e0 * k[16] * y[idx_H2] * y[idx_H]
                      + k[17] * y[idx_Hk] * y[idx_E]
                      - k[18] * y[idx_Hk] * y[idx_H] + 2.e0 * k[18] * y[idx_Hk] * y[idx_H]
                      - k[19] * y[idx_Hk] * y[idx_H] + 2.e0 * k[19] * y[idx_Hk] * y[idx_H]
                      + 2.e0 * k[20] * y[idx_Hk] * y[idx_Hj]
                      + 2.e0 * k[22] * y[idx_H2j] * y[idx_E]
                      + 2.e0 * k[23] * y[idx_H2j] * y[idx_E]
                      + k[24] * y[idx_H2j] * y[idx_Hk]
                      - 3.e0 * k[25] * y[idx_H] * y[idx_H] * y[idx_H] + k[25] * y[idx_H] * y[idx_H] * y[idx_H]
                      - 3.e0 * k[26] * y[idx_H] * y[idx_H] * y[idx_H] + k[26] * y[idx_H] * y[idx_H] * y[idx_H]
                      - 2.e0 * k[27] * y[idx_H2] * y[idx_H] * y[idx_H]
                      - 2.e0 * k[28] * y[idx_H2] * y[idx_H] * y[idx_H]
                      + k[29] * y[idx_Hj] * y[idx_D]
                      - k[30] * y[idx_H] * y[idx_Dj]
                      + k[33] * y[idx_H2] * y[idx_D]
                      + k[34] * y[idx_H2] * y[idx_D]
                      - k[35] * y[idx_HD] * y[idx_H],
                      - k[3] * y[idx_HE] * y[idx_E]
                      + k[4] * y[idx_HEj] * y[idx_E]
                      + k[5] * y[idx_HEj] * y[idx_E],
                      + k[9] * y[idx_Hk] * y[idx_H]
                      + k[10] * y[idx_Hk] * y[idx_H]
                      + k[13] * y[idx_H2j] * y[idx_H]
                      - k[14] * y[idx_H2] * y[idx_Hj]
                      - k[15] * y[idx_H2] * y[idx_E]
                      - k[16] * y[idx_H2] * y[idx_H]
                      + k[24] * y[idx_H2j] * y[idx_Hk]
                      + k[25] * y[idx_H] * y[idx_H] * y[idx_H]
                      + k[26] * y[idx_H] * y[idx_H] * y[idx_H]
                      - k[27] * y[idx_H2] * y[idx_H] * y[idx_H] + 2.e0 * k[27] * y[idx_H2] * y[idx_H] * y[idx_H]
                      - k[28] * y[idx_H2] * y[idx_H] * y[idx_H] + 2.e0 * k[28] * y[idx_H2] * y[idx_H] * y[idx_H]
                      - k[31] * y[idx_H2] * y[idx_Dj]
                      + k[32] * y[idx_HD] * y[idx_Hj]
                      - k[33] * y[idx_H2] * y[idx_D]
                      - k[34] * y[idx_H2] * y[idx_D]
                      + k[35] * y[idx_HD] * y[idx_H],
                      - k[29] * y[idx_Hj] * y[idx_D]
                      + k[30] * y[idx_H] * y[idx_Dj]
                      - k[33] * y[idx_H2] * y[idx_D]
                      - k[34] * y[idx_H2] * y[idx_D]
                      + k[35] * y[idx_HD] * y[idx_H]
                      - k[36] * y[idx_D] * y[idx_Hk]
                      + k[37] * y[idx_Dj] * y[idx_E],
                      + k[31] * y[idx_H2] * y[idx_Dj]
                      - k[32] * y[idx_HD] * y[idx_Hj]
                      + k[33] * y[idx_H2] * y[idx_D]
                      + k[34] * y[idx_H2] * y[idx_D]
                      - k[35] * y[idx_HD] * y[idx_H]
                      + k[36] * y[idx_D] * y[idx_Hk],
                      + k[0] * y[idx_H] * y[idx_E]
                      - k[1] * y[idx_Hj] * y[idx_E]
                      - k[2] * y[idx_Hj] * y[idx_E]
                      - k[11] * y[idx_H] * y[idx_Hj]
                      - k[12] * y[idx_H] * y[idx_Hj]
                      + k[13] * y[idx_H2j] * y[idx_H]
                      - k[14] * y[idx_H2] * y[idx_Hj]
                      - k[20] * y[idx_Hk] * y[idx_Hj]
                      - k[21] * y[idx_Hk] * y[idx_Hj]
                      - k[29] * y[idx_Hj] * y[idx_D]
                      + k[30] * y[idx_H] * y[idx_Dj]
                      + k[31] * y[idx_H2] * y[idx_Dj]
                      - k[32] * y[idx_HD] * y[idx_Hj],
                      + k[3] * y[idx_HE] * y[idx_E]
                      - k[4] * y[idx_HEj] * y[idx_E]
                      - k[5] * y[idx_HEj] * y[idx_E]
                      - k[6] * y[idx_HEj] * y[idx_E]
                      + k[7] * y[idx_HEjj] * y[idx_E],
                      + k[11] * y[idx_H] * y[idx_Hj]
                      + k[12] * y[idx_H] * y[idx_Hj]
                      - k[13] * y[idx_H2j] * y[idx_H]
                      + k[14] * y[idx_H2] * y[idx_Hj]
                      + k[21] * y[idx_Hk] * y[idx_Hj]
                      - k[22] * y[idx_H2j] * y[idx_E]
                      - k[23] * y[idx_H2j] * y[idx_E]
                      - k[24] * y[idx_H2j] * y[idx_Hk],
                      + k[29] * y[idx_Hj] * y[idx_D]
                      - k[30] * y[idx_H] * y[idx_Dj]
                      - k[31] * y[idx_H2] * y[idx_Dj]
                      + k[32] * y[idx_HD] * y[idx_Hj]
                      - k[37] * y[idx_Dj] * y[idx_E],
                      + k[6] * y[idx_HEj] * y[idx_E] - k[7] * y[idx_HEjj] * y[idx_E],
                      0.0, 0.0, 0.0, 0.0])

    # Arbitrary networks for usage in KROME
    elif (example=="KROME"):

        exec(open(kromefile).read())

    return y, rhs, nvector, neq, maxorder, nrea

# replaces the pragmas with fortran syntax


def ReplacePragmas(fh_list, fout):
    """
        Replaces pragmas in the *.F90 files. This includes writing out the
        different Jacobian or the LU-matrix in fortran syntax and writing
        down the fixed vector length.

    Parameters
    ----------
    fh_list : list of files
        List with two files which are read in. A common template file
        and a solver template file. These files provide pragmas and are
        read in.

    fout : list of files
        List with two files which are written out. A common template file
        and a solver template file. These files have substituted pragmas
        from the ingoing template files.

    """
    for k, fh in enumerate(fh_list):
        for row in fh:
            srow = row.strip()
            if(srow == "#ODEVEC_L"):
                if(np.shape(P)[0] < args.maxsize):
                    for i in range(L.shape[0]):
                        for j in range(L.shape[1]):
                            fout[k].write("      L(i," + str(i + 1) + "," + str(j + 1) + ") = " +
                                          fcode(L[i, j], source_format='free', standard=95) + "\n")
            elif(srow == "#ODEVEC_U"):
                if(np.shape(P)[0] < args.maxsize):
                    for i in range(U.shape[0]):
                        for j in range(U.shape[1]):
                            fout[k].write("      U(i," + str(i + 1) + "," + str(j + 1) + ") = " +
                                          fcode(U[i, j], source_format='free', standard=95) + "\n")
            elif(srow == "#ODEVEC_JAC"):
                if (args.packaging=="DENSE"):
                    for i in range(jac.shape[0]):
                        for j in range(jac.shape[1]):
                            fout[k].write("      jac(i," + str(i + 1) + "," + str(j + 1) + ") = " +
                                          fcode(P[i, j], source_format='free', standard=95) + "\n")
                elif (args.packaging=="COO"):
                    for i in range(nnz):
                        fout[k].write("     jac%ind1(i,"+ str(i+1) + ") = " +
                                          fcode(P.col_list()[i][0], source_format='free', standard=95) + "\n")
                        fout[k].write("     jac%ind2(i,"+ str(i+1) + ") = " +
                                          fcode(P.col_list()[i][1], source_format='free', standard=95) + "\n")
                        fout[k].write("     jac%values(i,"+ str(i+1) + ") = " +
                                          fcode(P.col_list()[i][2], source_format='free', standard=95) + "\n")
                elif (args.packaging=="CSC"):
                    for i in range(nnz):
                        fout[k].write("     jac%ind1(i,"+ str(i+1) + ")")
            elif(srow == "#ODEVEC_RHS"):
                for i in range(rhs.shape[0]):
                    fout[k].write("      rhs(i," + str(i + 1) + ") = " +
                                  fcode(rhs[i], source_format='free', standard=95) + "\n")
            elif(srow == "#ODEVEC_LU_PRESENT"):
                if(np.shape(P)[0] < args.maxsize):
                    fout[k].write("      logical :: LU_PRESENT = .TRUE.")
                else:
                    fout[k].write("      logical :: LU_PRESENT = .FALSE.")
            elif(srow == "#ODEVEC_VECTORLENGTH"):
                fout[k].write( "    integer :: nvector=" + str(nvector) + "\n")
            elif(srow == "#ODEVEC_COMMON_MODULE"):
                fout[k].write("      use" + "\n")
            elif(srow == "#ODEVEC_EQUATIONS"):
                fout[k].write("    integer :: neq=" + str(neq) + "\n")
            elif(srow == "#ODEVEC_MAXORDER"):
                fout[k].write("    integer :: maxorder=" + str(maxorder) + "\n")
            elif(srow == "#ODEVEC_NNZ"):
                fout[k].write("    integer :: nnz=" + str(nnz) + "\n")
            elif(srow == "#ODEVEC_NNZ_L"):
                fout[k].write("    integer :: nnz_l=" + str(nnz_l) + "\n")
            elif(srow == "#ODEVEC_NNZ_U"):
                fout[k].write("    integer :: nnz_u=" + str(nnz_u) + "\n")
            elif(srow == "#ODEVEC_REACTIONS"):
                fout[k].write("    integer :: nrea=" + str(nrea) + "\n")
            elif(srow == "#ODEVEC_DT_MIN"):
                fout[k].write( "    this%dt_min =" + str(args.dt_min) + "\n")
            else:
                srow = row.strip()
                fout[k].write(row)


def ReorderSystemCMK(P,y,rhs):
    """
    BUG: Reordering is not working properly at the moment.
    BUG: Nothing done with y_tmp, rhs_tmp.

    Reorders the Jacobian in order to have a better LU matrix
    with the CMK algorithm.

    Parameters
    ----------
    P : sympy matrix
        Matrix to be reordered in order to have a better.
    y : sympy matrix
        Function symbols of equation system.
    rhs : sympy matrix
        Right-hand-side of the problem.

    Returns
    -------
    P_order : sympy matrix
        Reordered sympy matrix of ingoing matrix P.
    y_tmp : sympy matrix
        Reordered sympy matrix of ingoing matrix y.
    rhs_tmp : sympy matrix
        Reorodered sympy matrix of ingoin matrix rhs.
    perm : scipy matrix
        Permuttation matrix.
    """

    # replace non-zero entries with ones for ordering algorithms in scipy
    P1 = np.array(P.col_list())
    P1[:, 2] = 1

    data = P1[:, 2]
    i = P1[:, 1]
    j = P1[:, 0]

    Psci = sparse.csc_matrix((data.astype(int), (i, j)))
    # apply Cuthill-McKee-Algorithm
    perm = sparse.csgraph.reverse_cuthill_mckee(Psci, symmetric_mode=False)

    P_order = zeros(Psci.shape[0],Psci.shape[1])
    y_tmp = y
    rhs_tmp = rhs
    for i in range(Psci.shape[0]):
        for j in range(Psci.shape[1]):
            P_order[i,j] = P[perm[i],perm[j]]

    print("\nFinished reordering with CMK-Algorithm...")
    print("Permutation list", perm, "\n")

    return P_order,y_tmp,rhs_tmp,perm


def ReorderSystemInvert(P,y,rhs):
    """
    BUG: Reordering is not working properly at the moment.

    Reorders the Jacobian in order to have a better LU matrix by
    inverting the matrix (graphically).
    """

    perm = range(len(rhs))[::-1]

    P_order = zeros(neq,neq)
    rhs_order = rhs
    y_order = y
    for i in range(neq):
        for j in range(neq):
            P_order[i,j] = P[perm[i],perm[j]]

    print("\nFinished reordering with simple inversed numbering...")
    print("Permutation list", perm, "\n")
    return [P_order,y_order,rhs_order,perm]


# command line execution
if __name__ == '__main__':
    """ Command line execution of preprocessor

    This is the standard procedure in order to run the preprocessor.
    The routine calls all necessary routines for precomputation, preordering,
    etc. It finally writes out the fortran files, which can be compiled in a
    next step.
    """
    commit = check_output("git rev-parse --short HEAD", shell=True).rstrip()
    print("#------------------------------------------------------------------#")
    print("#        Running OdeVec Preprocessing                              #")
    print("#        Commit Hash: " + str(commit) + "                                   #")
    print("#------------------------------------------------------------------#")

    # ---------------- parsing -----------------------------------------------#
    # Please note: ordering alphabatically with parsing argument
    # ------------------------------------------------------------------------#
    parser = argparse.ArgumentParser(
        description='Preprocessor for OdeVec. A vectorized ODE-solver built for high throughput.')
    parser.add_argument(
        '--checksparsity',
        default=False,
        help='Adds additional output for sparsity information. Attention: needs \
              quite alot of time, because the LU-decomposition is applied.')
    parser.add_argument(
        '--commonfile',
        default='src/odevec_commons.F90',
        help='path to the not preprocessed common file')
    parser.add_argument(
        '--commonout',
        default='build/odevec_commons.f90',
        help='path to the output common file')
    parser.add_argument(
        '--example',
        default="PRIMORDIAL",
        help='pre-defined networks to solve')
    parser.add_argument(
        '--krome_setupfile',
        default=None,
        help='path to an extra file from krome providing the ode')
    parser.add_argument(
        '--maxsize',
        type=int,
        default=5,
        help='maximum size of the Jacobian, which should still be symbolically LU-decomposed')
    parser.add_argument(
        '--dt_min',
        default=1e-10,
        required=False,
        help='change minimal timestep')
    parser.add_argument(
        '--nnz',
        default=0,
        help='Set the number of non-zeroes of the Jacobian. This option is generally not needed and set on its own')
    parser.add_argument(
        '--nnz_l',
        default=0,
        help='Set the number of non-zeroes of the L matrix. If not set the preprocessor will evaluate the nonzeroes on its own, which might take some time.')
    parser.add_argument(
        '--nnz_u',
        default=0,
        help='Set the number of non-zeroes of the U matrix. If not set the preprocessor will evaluate the nonzeroes on its own, which might take some time.')
    parser.add_argument(
        '--nvector',
        default=1,
        type=int,
        required=True,
        help='set vector length for the solver')
    parser.add_argument(
        '--ordering',
        default=None,
        help='preorder algorithm with a common algorithm')
    parser.add_argument(
        '--packaging',
        default="DENSE",
        help='Adds a sparse packaging format like CSC/CSR to the solver.')
    parser.add_argument(
        '--solverfile',
        default='src/odevec.F90',
        help='path to the not preprocessed solver file')
    parser.add_argument(
        '--solverout',
        default='build/odevec.f90',
        help='path to the solver output file')
    parser.add_argument(
        '--sparsity_structure',
        default=False,
        help='Print sparsity structure of the matrices to stdout.')

    args = parser.parse_args()

    if (args.krome_setupfile != None):
        args.example="KROME"

    # get right-hand-side and other
    y, rhs, nvector, neq, maxorder, nrea = \
            GetSystemToSolve(args.nvector,example=args.example,kromefile=args.krome_setupfile)

    # calculate jacobian
    jac = rhs.jacobian(y)

    # calculate P (used within BDF-methods)
    dt, beta = symbols('dt beta')
    P = SparseMatrix(eye(neq) - dt * beta * jac)

    # Reorder the system if chosen so
    if(args.ordering=="CMK"):
        [P, y, rhs, perm] = ReorderSystemCMK(P,y,rhs)
    if(args.ordering=="INVERT"):
        [P, y, rhs, perm] = ReorderSystemInvert(P,y,rhs)

    # only run symbolic LU-decompostion if the system is not too large
    maxsize = args.maxsize
    if(np.shape(P)[0] < maxsize or args.nnz_l == 0 or args.nnz_u == 0):
        print("#  Evaluating L and U matrices and check for sparsity structure...")
        #Piv, L, D, U = P.LUdecompositionFF()
        L,U,_ = P.LUdecomposition()

    # save information for sparsity structure
    if(args.nnz!=0):
        nnz = args.nnz
    else:
        nnz = P.nnz()
    if(args.nnz_l!=0):
        nnz_l = args.nnz_l
    else:
        nnz_l = L.nnz()
    if(args.nnz_u!=0):
        nnz_u = args.nnz_u
    else:
        nnz_u = U.nnz()

    print("#  NNZs in Jacobian: " + str(nnz) + " | Sparsity: " + str(float(nnz)/(neq*neq)))
    print("#  NNZs in L: " + str(nnz_l) + "        | Sparsity: " + str(float(nnz_l)/(neq*neq)))
    print("#  NNZs in U: " + str(nnz_u) + "        | Sparsity: " + str(float(nnz_u)/(neq*neq)))
    print("#  NNZs in LU: " + str(nnz_l+nnz_u) + "      | Sparsity: " + str(float(nnz_l+nnz_u)/(neq*neq)))

    # Choose packaging
    if(args.packaging=="CSC"):
        print("#  Packaging format: CSC (compressed sparse column)")
    elif(args.packaging=="COO"):
        print("#  Packaging format: COO (coordinate list packaging)")
    elif(args.packaging=="CSR"):
        print("#  Packaging format: CSR (compressed sparse row)")

    # write changes to file
    fh_list = [open(args.solverfile), open(args.commonfile)]
    fout = [open(args.solverout, "w"), open(args.commonout, "w")]

    # load example
    if(args.example=="ROBER"):
        copyfile("tests/rober.f90","build/test.f90")
    elif(args.example=="PRIMORDIAL"):
        copyfile("tests/primordial.f90","build/test.f90")
    elif(args.example=="KROME"):
        pass

    # search for pragmas and replace them with the correct
    print("#  Replacing pragmas in Fortran source code...")
    ReplacePragmas(fh_list, fout)
    for fout_single in fout: print("#  Wrote out file: " + fout_single.name)

    # some command line output
    if(args.sparsity_structure=="True"):
        print("#  Sparsity structures...")
        P.print_nonzero()
        L.print_nonzero()
        U.print_nonzero()
    if(np.shape(P)[0] < maxsize):
        print("sparsity structure of L:")
        L.print_nonzero()
    if(np.shape(P)[0] < maxsize):
        print("sparsity structure of U:")
        U.print_nonzero()

    print("#------------------------------------------------------------------#")
    print("#        OdeVec preprocessing done!                                #")
    print("#------------------------------------------------------------------#")
