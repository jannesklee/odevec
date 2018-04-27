#!/usr/bin/env python

from __future__ import print_function
from sympy import symbols, Matrix, eye, zeros, fcode, SparseMatrix, lambdify, nsimplify
from scipy import sparse
from shutil import copyfile
import numpy as np
import argparse

# sympy part, calculate L and U


def GetSystemToSolve(nvector,example):
    # Define here the system that you wish to solve. Below is Robertsons examples.
    # need to include here a function which gets the rhs by the inserting krome
    # krome file or directly by krome instead of by hand.
    # General indices
    nvector = nvector
    print("Used vector length nvector:", nvector)
    maxorder = 5

    #------------------------------------------------------------------------------#
    if (example=="ROBER"):
        # robertson's test
        neq = 3
        y = symbols('y(i\,1:%d)'%(neq+2))
        k = 0
        nrea = 0

        rhs = Matrix([-0.04*y[0]+1e4*y[1]*y[2],0.04*y[0]-3e7*y[1]*y[1]-1e4*y[1]*y[2],3e7*y[1]*y[1]])

    elif (example=="PRIMORDIAL"):
        # Primordial network with 16 different species
        neq = 16
        nrea = 38


        idx_E = 1
        idx_Hk = 2
        idx_H = 3
        idx_HE = 4
        idx_H2 = 5
        idx_D = 6
        idx_HD = 7
        idx_Hj = 8
        idx_HEj = 9
        idx_H2j = 10
        idx_Dj = 11
        idx_HEjj = 12
        idx_CR = 13
        idx_g = 14
        idx_Tgas = 15
        idx_dummy = 16

        y = symbols('y(i\,0:%d)' % (neq + 1))
        k = symbols('k(i\,0:%d)' % (nrea + 1))

        # RHS of primordial network
        rhs = Matrix([- k[1] * y[idx_H] * y[idx_E] + 2.e0 * k[1] * y[idx_H] * y[idx_E] - k[2] * y[idx_Hj] * y[idx_E] - k[3] * y[idx_Hj] * y[idx_E] - k[4] * y[idx_HE] * y[idx_E] + 2.e0 * k[4] * y[idx_HE] * y[idx_E] - k[5] * y[idx_HEj] * y[idx_E] - k[6] * y[idx_HEj] * y[idx_E] - k[7] * y[idx_HEj] * y[idx_E] + 2.e0 * k[7] * y[idx_HEj] * y[idx_E] - k[8] * y[idx_HEjj] * y[idx_E] - k[9] * y[idx_H] * y[idx_E] + k[10] * y[idx_Hk] * y[idx_H] + k[11] * y[idx_Hk] * y[idx_H] - k[16] * y[idx_H2] * y[idx_E] + k[16] * y[idx_H2] * y[idx_E] - k[18] * y[idx_Hk] * y[idx_E] + 2.e0 * k[18] * y[idx_Hk] * y[idx_E] + k[19] * y[idx_Hk] * y[idx_H] + k[20] * y[idx_Hk] * y[idx_H] + k[22] * y[idx_Hk] * y[idx_Hj] - k[23] * y[idx_H2j] * y[idx_E] - k[24] * y[idx_H2j] * y[idx_E] + k[37] * y[idx_D] * y[idx_Hk] - k[38] * y[idx_Dj] * y[idx_E], + k[9] * y[idx_H] * y[idx_E] - k[10] * y[idx_Hk] * y[idx_H] - k[11] * y[idx_Hk] * y[idx_H] - k[18] * y[idx_Hk] * y[idx_E] - k[19] * y[idx_Hk] * y[idx_H] - k[20] * y[idx_Hk] * y[idx_H] - k[21] * y[idx_Hk] * y[idx_Hj] - k[22] * y[idx_Hk] * y[idx_Hj] - k[25] * y[idx_H2j] * y[idx_Hk] - k[37] * y[idx_D] * y[idx_Hk], - k[1] * y[idx_H] * y[idx_E] + k[2] * y[idx_Hj] * y[idx_E] + k[3] * y[idx_Hj] * y[idx_E] - k[9] * y[idx_H] * y[idx_E] - k[10] * y[idx_Hk] * y[idx_H] - k[11] * y[idx_Hk] * y[idx_H] - k[12] * y[idx_H] * y[idx_Hj] - k[13] * y[idx_H] * y[idx_Hj] - k[14] * y[idx_H2j] * y[idx_H] + k[15] * y[idx_H2] * y[idx_Hj] + 2.e0 * k[16] * y[idx_H2] * y[idx_E] - k[17] * y[idx_H2] * y[idx_H] + 3.e0 * k[17] * y[idx_H2] * y[idx_H] + k[18] * y[idx_Hk] * y[idx_E] - k[19] * y[idx_Hk] * y[idx_H] + 2.e0 * k[19] * y[idx_Hk] * y[idx_H] - k[20] * y[idx_Hk] * y[idx_H] + 2.e0 * k[20] * y[idx_Hk] * y[idx_H] + 2.e0 * k[21] * y[idx_Hk] * y[idx_Hj] + 2.e0 * k[23] * y[idx_H2j] * y[idx_E] + 2.e0 * k[24] * y[idx_H2j] * y[idx_E] + k[25] * y[idx_H2j] * y[idx_Hk] - 3.e0 * k[26] * y[idx_H] * y[idx_H] * y[idx_H] + k[26] * y[idx_H] * y[idx_H] * y[idx_H] - 3.e0 * k[27] * y[idx_H] * y[idx_H] * y[idx_H] + k[27] * y[idx_H] * y[idx_H] * y[idx_H] - 2.e0 * k[28] * y[idx_H2] * y[idx_H] * y[idx_H] - 2.e0 * k[29] * y[idx_H2] * y[idx_H] * y[idx_H] + k[30] * y[idx_Hj] * y[idx_D] - k[31] * y[idx_H] * y[idx_Dj] + k[34] * y[idx_H2] * y[idx_D] + k[35] * y[idx_H2] * y[idx_D] - k[36] * y[idx_HD] * y[idx_H], - k[4] * y[idx_HE] * y[idx_E] + k[5] * y[idx_HEj] * y[idx_E] + k[6] * y[idx_HEj] * y[idx_E], + k[10] * y[idx_Hk] * y[idx_H] + k[11] * y[idx_Hk] * y[idx_H] + k[14] * y[idx_H2j] * y[idx_H] - k[15] * y[idx_H2] * y[idx_Hj] - k[16] * y[idx_H2] * y[idx_E] - k[17] * y[idx_H2] * y[idx_H] + k[25] * y[idx_H2j] * y[idx_Hk] + k[26] * y[idx_H] * y[idx_H] * y[idx_H] + k[27] * y[idx_H] * y[idx_H] * y[idx_H] - k[28] * y[idx_H2] * y[idx_H] * y[idx_H] + 2.e0 * k[28] * y[idx_H2] * y[idx_H] * y[idx_H] - k[29] * y[idx_H2] * y[idx_H] * y[idx_H] + 2.e0 * k[29] * y[idx_H2] * y[idx_H] * y[idx_H] - k[32] * y[idx_H2] * y[idx_Dj] + k[33] * y[idx_HD] * y[idx_Hj] - k[34] * y[idx_H2] * y[idx_D] - k[35] * y[idx_H2] * y[idx_D] + k[36] * y[idx_HD] * y[idx_H], - k[30] * y[idx_Hj] * y[idx_D] + k[31] * y[idx_H] * y[idx_Dj] - k[34] * y[idx_H2] * y[idx_D] - k[35] * y[idx_H2] * y[idx_D] + k[36] * y[idx_HD] * y[idx_H] - k[37] * y[idx_D] * y[idx_Hk] + k[38] * y[idx_Dj] * y[idx_E], + k[32] * y[idx_H2] * y[idx_Dj] - k[33] * y[idx_HD] * y[idx_Hj] + k[34] * y[idx_H2] * y[idx_D] + k[35] * y[idx_H2] * y[idx_D] - k[36] * y[idx_HD] * y[idx_H] + k[37] * y[idx_D] * y[idx_Hk], + k[1] * y[idx_H] * y[idx_E] - k[2] * y[idx_Hj] * y[idx_E] - k[3] * y[idx_Hj] * y[idx_E] - k[12] * y[idx_H] * y[idx_Hj] - k[13] * y[idx_H] * y[idx_Hj] + k[14] * y[idx_H2j] * y[idx_H] - k[15] * y[idx_H2] * y[idx_Hj] - k[21] * y[idx_Hk] * y[idx_Hj] - k[22] * y[idx_Hk] * y[idx_Hj] - k[30] * y[idx_Hj] * y[idx_D] + k[31] * y[idx_H] * y[idx_Dj] + k[32] * y[idx_H2] * y[idx_Dj] - k[33] * y[idx_HD] * y[idx_Hj], + k[4] * y[idx_HE] * y[idx_E] - k[5] * y[idx_HEj] * y[idx_E] - k[6] * y[idx_HEj] * y[idx_E] - k[7] * y[idx_HEj] * y[idx_E] + k[8] * y[idx_HEjj] * y[idx_E], + k[12] * y[idx_H] * y[idx_Hj] + k[13] * y[idx_H] * y[idx_Hj] - k[14] * y[idx_H2j] * y[idx_H] + k[15] * y[idx_H2] * y[idx_Hj] + k[22] * y[idx_Hk] * y[idx_Hj] - k[23] * y[idx_H2j] * y[idx_E] - k[24] * y[idx_H2j] * y[idx_E] - k[25] * y[idx_H2j] * y[idx_Hk], + k[30] * y[idx_Hj] * y[idx_D] - k[31] * y[idx_H] * y[idx_Dj] - k[32] * y[idx_H2] * y[idx_Dj] + k[33] * y[idx_HD] * y[idx_Hj] - k[38] * y[idx_Dj] * y[idx_E], + k[7] * y[idx_HEj] * y[idx_E] - k[8] * y[idx_HEjj] * y[idx_E], 0.0, 0.0, 0.0, 0.0])
        #------------------------------------------------------------------------------#

    return [y, rhs, nvector, neq, maxorder, nrea]

# replaces the pragmas with fortran syntax


def ReplacePragmas(fh_list, fout):
    for k, fh in enumerate(fh_list):
        for row in fh:
            srow = row.strip()
            if(srow == "#ODEVEC_L"):
                if(np.shape(P)[0] < 5):
                    for i in range(L.shape[0]):
                        for j in range(L.shape[1]):
                            fout[k].write("      L(i," + str(i + 1) + "," + str(j + 1) + ") = " +
                                          fcode(L[i, j], source_format='free', standard=95) + "\n")
            elif(srow == "#ODEVEC_U"):
                if(np.shape(P)[0] < 5):
                    for i in range(U.shape[0]):
                        for j in range(U.shape[1]):
                            fout[k].write("      U(i," + str(i + 1) + "," + str(j + 1) + ") = " +
                                          fcode(U[i, j], source_format='free', standard=95) + "\n")
            elif(srow == "#ODEVEC_JAC"):
                for i in range(jac.shape[0]):
                    for j in range(jac.shape[1]):
                        fout[k].write("      jac(i," + str(i + 1) + "," + str(j + 1) + ") = " +
                                      fcode(P[i, j], source_format='free', standard=95) + "\n")
            elif(srow == "#ODEVEC_RHS"):
                for i in range(rhs.shape[0]):
                    fout[k].write("      rhs(i," + str(i + 1) + ") = " +
                                  fcode(rhs[i], source_format='free', standard=95) + "\n")
            elif(srow == "#ODEVEC_LU_PRESENT"):
                if(np.shape(P)[0] < 5):
                    fout[k].write("      LOGICAL :: LU_PRESENT = .TRUE.")
                else:
                    fout[k].write("      LOGICAL :: LU_PRESENT = .FALSE.")
            elif(srow == "#ODEVEC_VECTORLENGTH"):
                fout[k].write(
                    "      integer :: nvector=" +
                    str(nvector) +
                    "\n")
            elif(srow == "#ODEVEC_COMMON_MODULE"):
                fout[k].write("      use" + "\n")
            elif(srow == "#ODEVEC_EQUATIONS"):
                fout[k].write("      integer :: neq=" + str(neq) + "\n")
            elif(srow == "#ODEVEC_MAXORDER"):
                fout[k].write(
                    "      integer :: maxorder=" +
                    str(maxorder) +
                    "\n")
            elif(srow == "#ODEVEC_REACTIONS"):
                fout[k].write("      integer :: nrea=" + str(nrea) + "\n")
            else:
                srow = row.strip()
                fout[k].write(row)


def ReorderSystemCMK(P,y,rhs):
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

    return [P_order,y_tmp,rhs_tmp,perm]

def ReorderSystemInvert(P,y,rhs):

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
    print("#----------------------------------------------#")
    print("Running OdeVec Preprocessing")
    print("#----------------------------------------------#")
    # parsing
    parser = argparse.ArgumentParser(
        description='Preprocessor for OdeVec. A vectorized ODE-solver built for high throughput.')
    parser.add_argument(
        '--solverfile',
        default='src/odevec.F90',
        help='path to the not preprocessed solver file')
    parser.add_argument(
        '--solverout',
        default='build/odevec.f90',
        help='path to the solver output file')
    parser.add_argument(
        '--commonfile',
        default='src/odevec_commons.F90',
        help='path to the not preprocessed common file')
    parser.add_argument(
        '--commonout',
        default='build/odevec_commons.f90',
        help='path to the output common file')
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
        '--example',
        default="ROBER",
        help='pre-defined networks to solve')
    parser.add_argument(
        '--maxsize',
        default=5,
        help='maximum size of the Jacobian, which should still be symbolically LU-decomposed')
    args = parser.parse_args()

    # get right-hand-side
    [y, rhs, nvector, neq, maxorder, nrea] = GetSystemToSolve(args.nvector,example=args.example)

    # calculate jacobian
    jac = nsimplify(SparseMatrix(rhs.jacobian(y)[:, 1:]))
    # calculate P (used within BDF-methods)
    dt, beta = symbols('dt beta')
    P = SparseMatrix(eye(neq) - dt * beta * jac)

    # Reorder the system if chosen so
    if(args.ordering=="CMK"):
        [P, y, rhs, perm] = ReorderSystemCMK(P,y,rhs)
    if(args.ordering=="INVERT"):
        [P, y, rhs, perm] = ReorderSystemInvert(P,y,rhs)

    maxsize = args.maxsize
    # only run symbolic LU-decompostion if the system is not too large
    if(np.shape(P)[0] < maxsize):
        #L,U,_ = P.LUdecomposition()
        Piv, L, D, U = P.LUdecompositionFF()
#        print(U.nnz())
#        print(L.RL)

    # write changes to file
    fh_list = [open(args.solverfile), open(args.commonfile)]
    fout = [open(args.solverout, "w"), open(args.commonout, "w")]

    if(args.example=="ROBER"):
        copyfile("tests/rober.f90","build/test.f90")
    elif(args.example=="PRIMORDIAL"):
        copyfile("tests/primordial.f90","build/test.f90")

    # search for pragmas and replace them with the correct
    ReplacePragmas(fh_list, fout)

    # some command line output
    print("Sparsity structure of P:")
    P.print_nonzero()
    if(np.shape(P)[0] < maxsize):
        print("sparsity structure of L:")
        L.print_nonzero()
    if(np.shape(P)[0] < maxsize):
        print("sparsity structure of U:")
        U.print_nonzero()

    print("OdeVec preprocessing done!")
