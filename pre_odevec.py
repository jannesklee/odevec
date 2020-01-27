#!/usr/bin/env python
###############################################################################
#                                                                             #
#    Copyright (C) 2019                                                       #
#    Jannes Klee <jklee@astrophysik.uni-kiel.de>                              #
#                                                                             #
#    This file is part of OdeVec.                                             #
#                                                                             #
#    OdeVec is free software: you can redistribute it and/or modify           #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    OdeVec is distributed in the hope that it will be useful,                #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with OdeVec.  If not, see <https://www.gnu.org/licenses/>.         #
#                                                                             #
###############################################################################

""" Python preprocessor for the stiff BDF solver OdeVec

This script sets up the system to solve and does some precomputations
if possible. This includes writing down Fortran code in *.F90 template
files. There, placeholder (Pragmas) are placed, which are searched and
replaced. The precomputations includes the symbolic computation of the
Jacobian and preordering of the latter.
"""

from __future__ import print_function
from sympy import *
from scipy import sparse
from shutil import copyfile
from subprocess import check_output
import sys
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
    print("#  Setting up system..")
    print("#  Vector length \"nvector\":", nvector)

    #------------------------------------------------------------------------------#
    # Robertson's example
    if (example=="ROBER"):
        print("#  Setup: ROBER test")

        neq = 3
        k = 0
        nrea = 0

        # define sympy symbols
        y = list(symbols('y(\:\,1:%d)'%(neq+1)))

        rhs = Matrix([-0.04*y[0]+1e4*y[1]*y[2],
                      0.04*y[0]-3e7*y[1]*y[1]-1e4*y[1]*y[2],
                      3e7*y[1]*y[1]])

    # Primordial network with 16 different species
    elif (example=="PRIMORDIAL"):
        print("#  Setup: PRIMORDIAL test")

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
        y = list(symbols('y(\:\,1:%d)' % (neq + 1)))
        k = list(symbols('k(\:\,1:%d)' % (nrea + 1)))

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
    elif (example=="OREGO"):
        print("#  Setup: OREGO test")

        neq = 3
        nrea = 0

        # define sympy symbols
        y = list(symbols('y(\:\,1:%d)'%(neq+1)))

        k1 = 77.27
        k2 = 8.375e-6
        k3 = 0.161

        rhs = Matrix([k1*y[0]+k1*y[1]-k1*k2*y[0]**2-k1*y[0]*y[1],
                      1.0/k1*y[2]-1.0/k1*y[1]-1.0/k1*y[0]*y[1],
                      k3*y[0]-k3*y[2]])

    # Arbitrary networks for usage in KROME
    elif (example=="KROME"):
        print("#  Setup: KROME file read in")

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
            if(srow == "#ODEVEC_LU"):
                if(args.LUmethod==2):
                    if (args.packaging=="DENSE"):
                        for i in range(LU.shape[0]):
                            for j in range(LU.shape[1]):
                                if (not LU[i,j] == 0.0):
                                    fout[k].write("      LU(:," + str(i + 1) + "," + str(j + 1) + ") = " +
                                                  fcode(LU[i, j], source_format='free', standard=95) + "\n")
                    elif (args.packaging=="CSC"):
                        for i in range(LU.nnz()):
                            fout[k].write("      LU%sdata(:," + str(i + 1) + ") = " +
                                          fcode(LU.col_list()[i][2], source_format='free', standard=95) + "\n")
            elif(srow == "#ODEVEC_JAC"):
                if (args.packaging=="DENSE"):
                    for j in range(jac.shape[1]):
                        for i in range(jac.shape[0]):
                            if(args.heatcool==1):
                                if((i==neq-2 and j==neq-2)): #Tgas
                                        fout[k].write("      dyy(:) = y(:," + str(j+1) + ")*1d-2\n")
                                        fout[k].write("      yy(:," + str(j+1) + ") = y(:," + str(j+1) +") + dyy(:)\n")
                                        fout[k].write("      call GetRHS(this,yy(:,:),dy(:,:))\n")
                                        fout[k].write("      do i=1,this%neq\n")
                                        fout[k].write("        jac(:,i,"+ str(j+1) + ") = -beta*dt*(dy(:,i)/dyy(:))\n")
                                        fout[k].write("      end do\n")
                                        fout[k].write("      jac(:,"+ str(i+1) +","+ str(j+1) + ") = 1.0+jac(:,"+ str(i+1) +","+ str(j+1) + ")\n")
                                elif((i==neq-2) and (j!=neq-2) and (j<neq-4)): #Tgas
                                        fout[k].write("      yy(:,:) = y(:,:)\n")
                                        fout[k].write("      dyy(:) = y(:," + str(j+1) +")*1d-2\n")
                                        fout[k].write("      yy(:," + str(j+1) + ") = y(:," + str(j+1) +") + dyy(:)\n")
                                        fout[k].write("      dy1(:) = " + krome_heatcool_string + "\n")
                                        fout[k].write("      where (dyy>0d0)\n")
                                        fout[k].write("          jac(:,"+ str(i+1) + "," + str(j+1) + ") = -beta*dt*(dy1(:)-dy0(:))/dyy(:)\n")
                                        fout[k].write("      end where\n")
                                elif(not (i==neq-2) and (j==neq-2) ): #Tgas
                                        pass
                                else:
                                        fout[k].write("      jac(:," + str(i + 1) + "," + str(j + 1) + ") = " +
                                                      fcode(P_order[i, j], source_format='free', standard=95) + "\n")
                            else:
                                fout[k].write("      jac(:," + str(i + 1) + "," + str(j + 1) + ") = " +
                                              fcode(P_order[i, j], source_format='free', standard=95) + "\n")
                elif (args.packaging=="CSC"):
                    for i in range(P_order.nnz()):
                        if (P_order.col_list()[i][2] == fills):
                            fout[k].write("      jac%sdata(:," + str(i + 1) + ") = " +
                                          fcode(0.0, source_format='free', standard=95) + "\n")
                        else:
                            fout[k].write("      jac%sdata(:," + str(i + 1) + ") = " +
                                          fcode(P_order.col_list()[i][2], source_format='free', standard=95) + "\n")
            elif(srow == "#ODEVEC_RHS"):
                for i in range(rhs.shape[0]):
                    if((i==neq-2) and (args.heatcool==1)):
                        fout[k].write("      rhs(:," + str(i + 1) + ") = " + krome_heatcool_string0 + "\n")
                    else:
                        fout[k].write("      rhs(:," + str(i + 1) + ") = " +
                                      fcode(rhs[i], source_format='free', standard=95) + "\n")
            elif(srow == "#ODEVEC_DY0"):
                if(args.heatcool==1):
                        fout[k].write("      dy0(:) = " + krome_heatcool_string0 + "\n")
            elif(srow == "#ODEVEC_PERMUTATIONS"):
                fout[k].write( "    this%Perm = " +
                              fcode(Perm+1, source_format='free', standard=95) + "\n")
            elif(srow == "#ODEVEC_LU_METHOD"):
                fout[k].write("    INTEGER :: LUmethod = " + str(int(LUmethod)) + "\n")
            elif(srow == "#ODEVEC_ALLOCATE_LU"):
                if (args.packaging=="DENSE"):
                    fout[k].write( "    allocate(this%LU(this%nvector,this%neq,this%neq),STAT=err)\n" +
                                   "    this%LU(:,:,:) = 0.0")
                if (args.packaging=="CSC"):
                    fout[k].write( "    allocate( &\n\
              this%LU%sdata(this%nvector,this%nnz), &\n\
              this%LU%u_col_start(this%neq+1), &\n\
              this%LU%l_col_start(this%neq), &\n\
              this%LU%row_index(this%nnz), &\n\
              STAT=err)\n")
                    fout[k].write( "    this%LU%sdata(:,:) = 0.0\n")
            elif(srow == "#ODEVEC_DEALLOCATE_LU"):
                if (args.packaging=="DENSE"):
                    fout[k].write( "    deallocate(this%LU) \n")
                if (args.packaging=="CSC"):
                    fout[k].write( "    deallocate( & \n \
              this%LU%sdata, &\n \
              this%LU%u_col_start, &\n \
              this%LU%l_col_start, &\n \
              this%LU%row_index)\n")
            elif(srow == "#ODEVEC_SET_LU_SPARSITY"):
                if(args.packaging=="CSC"):
                    fout[k].write( "    this%LU%row_index = " +
                            fcode(np.array(row_index[:])+1, source_format='free', standard=95) + "\n")
                    fout[k].write( "    this%LU%u_col_start = " +
                            fcode(np.array(u_col_start[:])+1, source_format='free', standard=95) + "\n")
                    fout[k].write( "    this%LU%l_col_start = " +
                            fcode(np.array(l_col_start[:])+1, source_format='free', standard=95) + "\n")
            elif(srow == "#ODEVEC_VECTORLENGTH"):
                fout[k].write( "    integer :: nvector=" + str(nvector) + "\n")
            elif(srow == "#ODEVEC_LU_MATRIX"):
                if (args.packaging=="DENSE"):
                    fout[k].write( "    double precision, pointer, " +
                                   "dimension(:,:,:) :: LU       !> jacobian & LU Matrix\n")
                elif (args.packaging=="CSC"):
                    fout[k].write( "    type(csc_matrix) :: LU      !> jacobian & LU Matrix\n")
            elif(srow == "#ODEVEC_COMMON_MODULE"):
                fout[k].write("      use" + "\n")
            elif(srow == "#ODEVEC_EQUATIONS"):
                fout[k].write("    integer :: neq=" + str(neq) + "\n")
            elif(srow == "#ODEVEC_MAXORDER"):
                fout[k].write("    integer :: maxorder=" + str(maxorder) + "\n")
            elif(srow == "#ODEVEC_NNZ"):
                fout[k].write("    integer :: nnz=" + str(LU.nnz()) + "\n")
            elif(srow == "#ODEVEC_PACKAGING"):
                if (args.packaging=="DENSE"):
                    fout[k].write("    character(len=10) :: packaging=\"dense\"\n")
                elif (args.packaging=="CSC"):
                    fout[k].write("    character(len=10) :: packaging=\"csc\"\n")
            elif(srow == "#ODEVEC_JAC_DEC"):
                if (args.packaging=="DENSE"):
                    fout[k].write("    double precision, dimension(this%nvector,this%neq,this%neq) :: jac\n")
                if (args.packaging=="CSC"):
                    fout[k].write("    type(csc_matrix) :: jac\n")
            elif(srow == "#ODEVEC_LU_DEC"):
                if (args.packaging=="DENSE"):
                    fout[k].write("    double precision, dimension(this%nvector,this%neq,this%neq) :: LU\n")
                if (args.packaging=="CSC"):
                    fout[k].write("    type(csc_matrix) :: LU\n")
            elif(srow == "#ODEVEC_REACTIONS"):
                fout[k].write("    integer :: nrea=" + str(nrea) + "\n")
            elif(srow == "#ODEVEC_DT_MIN"):
                fout[k].write( "   this%MinimalTimestep =" + str(args.dt_min) + "\n")
            else:
                srow = row.strip()
                fout[k].write(row)


def ReorderSystemCMK(P,y,rhs):
    """
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
    Perm : scipy matrix
        Permutation matrix.
    """

    # replace non-zero entries with ones for ordering algorithms in scipy
    P1 = np.array(P.col_list())
    P1[:, 2] = 1

    data = P1[:, 2]
    i = P1[:, 1]
    j = P1[:, 0]

    Psci = sparse.csc_matrix((data.astype(int), (i, j)))

    # apply Cuthill-McKee-Algorithm
    Perm = sparse.csgraph.reverse_cuthill_mckee(Psci, symmetric_mode=False)

    P_order = SparseMatrix.zeros(neq)
    rhs_order = rhs
    y_order = y
    for i in range(Psci.shape[0]):
        for j in range(Psci.shape[1]):
            P_order[i,j] = P[Perm[i],Perm[j]]

    print("#  Reordering with CMK-Algorithm..")
    print("#  Permutation list", Perm)

    return P_order,y_order,rhs_order,Perm


def ReorderSystemInvert(P,y,rhs):
    """
    Reorders the Jacobian in order to have a better LU matrix by
    inverting the matrix (graphically).
    """

    Perm = np.arange(len(rhs))[::-1]

    P_order = SparseMatrix.zeros(neq)
    rhs_order = rhs
    y_order = y
    for i in range(neq):
        for j in range(neq):
            P_order[i,j] = P[Perm[i],Perm[j]]

    print("#  Reordering: Inversed indices..")
    print("#  Permutation list", Perm)
    return P_order,y_order,rhs_order,Perm


def ReorderSystemFewestFirst(P,y,rhs):
    """
    Reorders the Jacobian by ordering after the combined amount
    of non-zeros along a column and rows.
    """

    Perm = np.arange(len(rhs))

    sorter = np.zeros(neq)
    for i in range(neq):
        sorter[i] = P[:,i].nnz() + P[i,:].nnz() - 1

    sort_indices = sorter.argsort()
    Perm = Perm[sort_indices]

    P_order = SparseMatrix.zeros(neq)
    rhs_order = rhs
    y_order = y
    for i in range(neq):
        for j in range(neq):
            P_order[i,j] = P[Perm[i],Perm[j]]

    print("#  Reordering: putting fewest reactions first..")

    return P_order, y_order, rhs_order, Perm

# command line execution
if __name__ == '__main__':
    """ Command line execution of preprocessor

    This is the standard procedure in order to run the preprocessor.
    The routine calls all necessary routines for precomputation, preordering,
    etc. It finally writes out the fortran files, which can be compiled in a
    next step.
    """
    commit = check_output("git rev-parse --short HEAD", shell=True).rstrip()
    dirty = check_output("git diff --quiet || echo 'dirty'", shell=True).rstrip()
    print("#------------------------------------------------------------------#")
    print("#        Running OdeVec Preprocessing                              #")
    print("#        Commit Hash: " + commit.decode("utf-8") +" - "+ \
            dirty.decode("utf-8") + "                              #")
    print("#------------------------------------------------------------------#")

    # ---------------- parsing -----------------------------------------------#
    # Please note: ordering alphabatically with parsing argument
    # ------------------------------------------------------------------------#
    parser = argparse.ArgumentParser(
            description='Preprocessor for OdeVec. A vectorized ODE-solver built for high throughput.')
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
            '--LUmethod',
            type=int,
            default=1,
            help='states how to calculate LU (Jacobian) method')
    parser.add_argument(
            '--dt_min',
            default=1e-10,
            required=False,
            help='change minimal timestep')
    parser.add_argument(
            '--nvector',
            default=1,
            type=int,
            required=True,
            help='set vector length for the solver')
    parser.add_argument(
            '--heatcool',
            default=0,
            type=int,
            required=False,
            help='set to 1 if cooling/heating should be solved within ODE in KROMEian way')
    parser.add_argument(
            '--equilibrium_H',
            default=0,
            type=int,
            required=False,
            help='set 1 in order to avoid problems for high number densities of hydrogen')
    parser.add_argument(
            '--ordering',
            default="None",
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

    if ((args.equilibrium_H==1) and args.packaging != "DENSE"):
        sys.exit("args.equilibrium_H is only allowed with DENSE packaging.")
    if ((args.equilibrium_H==1) and (args.ordering != "None")):
        sys.exit("args.equilibrium_H is only allowed without ordering.")
    if ((args.equilibrium_H==1) and args.LUmethod != 3):
        sys.exit("args.equilibrium_H is currently only usable with LUmethod==3.")

    # get right-hand-side and other
    y, rhs, nvector, neq, maxorder, nrea = \
            GetSystemToSolve(args.nvector,example=args.example,kromefile=args.krome_setupfile)

    # calculate jacobian
    jac = rhs.jacobian(y)

    # overwrite jacobian and rhs when cooling/heating is present in krome
    # \todo find a better solution for this ugly one
#    jac[idx_Tgas] = (heating(y,) - cooling(y)) * (krome_gamma - 1.0) / boltzmann_erg / Add(*[y[m] for m in range("+str(neq-4)+")])

    # calculate P (used within BDF-methods)
    dt, beta = symbols('dt beta')
    P = SparseMatrix(eye(neq) - dt * beta * jac)

    # Reorder the system if chosen so
    if(args.ordering=="CMK"):
        P_order, y_order, rhs_order, Perm = ReorderSystemCMK(P,y,rhs)
    elif(args.ordering=="INVERT"):
        P_order, y_order, rhs_order, Perm = ReorderSystemInvert(P,y,rhs)
    elif(args.ordering=="FF"):
        P_order, y_order, rhs_order, Perm = ReorderSystemFewestFirst(P,y,rhs)
    else:
        P_order = P
        y_order = y
        rhs_order = rhs
        Perm = np.arange(len(rhs))

    # only run symbolic LU-decompostion if the system is not too large
    LUmethod = args.LUmethod

    print("#  Evaluating L and U matrices and check for sparsity structure..")
    LU, Piv = P_order.LUdecomposition_Simple()

    # save information for sparsity structure
    print("#  NNZs in Jacobian: " + str(P.nnz()) + " | Sparsity: " + str(float(P.nnz())/(neq*neq)))
    print("#  NNZs in LU-Matrix: " + str(LU.nnz()) + " | Sparsity: " + str(float(LU.nnz())/(neq*neq)))

    # Choose packaging
    if(args.packaging=="CSC"):
        print("#  Packaging format: CSC (compressed sparse column)")
        # pad fillins to P matrix in order to have the same structure as LU
        P_order = Matrix(P_order)
        LU = Matrix(LU)
        fills = sympify("fills")                    # fill in placeholder
        for j in range(neq):
            for i in range(neq):
                if (P_order[i,j] == 0 and LU[i,j] != 0):
                    P_order[i,j] = fills
        P_order = SparseMatrix(P_order)
        LU = SparseMatrix(LU)
    elif(args.packaging=="DENSE"):
        print("#  Packaging format: dense (no compression)")
    elif(args.packaging=="COO"):
        sys.exit("# Error: COO-format (coordinate list packaging) NOT SUPPORTED")
    elif(args.packaging=="CSR"):
        sys.exit("# Error: CSR-format (compressed sparse row) NOT SUPPORTED")
    else:
        sys.exit("# Error: This format is not known.")


    # define the sparsity structure of LU-Matrix
    row_index = zeros(LU.nnz(),1)
    col_start = zeros(neq,1)
    u_col_start = zeros(neq+1,1)
    l_col_start = zeros(neq,1)
    value = zeros(LU.nnz(),1)

    j_old = 0
    for k in range(LU.nnz()):
        i = LU.col_list()[k][0]              # row index
        j = LU.col_list()[k][1]              # column index
        if(j-1==j_old):
            u_col_start[j] = k
        if(i==j):
            l_col_start[j] = k
        row_index[k] = i                     # row index
        j_old = j
    u_col_start[j+1] = LU.nnz()

    # write changes to file
    fh_list = [open(args.solverfile), open(args.commonfile)]
    fout = [open(args.solverout, "w"), open(args.commonout, "w")]

    # load example
    if(args.example=="ROBER"):
        copyfile("tests/rober.f90","build/test.f90")
    elif(args.example=="PRIMORDIAL"):
        copyfile("tests/primordial.f90","build/test.f90")
    elif(args.example=="OREGO"):
        copyfile("tests/orego.f90","build/test.f90")
    elif(args.example=="KROME"):
        pass

    krome_heatcool_string0 = "(heating(y(:,:), Tgas(:), k(:,:), nH2dust(:)) &\n \
             - cooling(y(:,:), Tgas(:)) #KROME_coolingQuench #KROME_coolfloor) &\n \
             * (krome_gamma(:) - 1.d0) / boltzmann_erg / sum(y(:,1:nmols),dim=2)"
    krome_heatcool_string = "(heating(yy(:,:), Tgas(:), k(:,:), nH2dust(:)) &\n \
             - cooling(yy(:,:), Tgas(:)) #KROME_coolingQuench #KROME_coolfloor) &\n \
             * (krome_gamma(:) - 1.d0) / boltzmann_erg / sum(y(:,1:nmols),dim=2)"

    # search for pragmas and replace them with the correct
    print("#  Replacing pragmas in Fortran source code..")
    ReplacePragmas(fh_list, fout)
    for fout_single in fout: print("#  Wrote out file: " + fout_single.name)

    # some command line output
    if(args.sparsity_structure=="True"):
        print("#  Sparsity structures..")
        print("#  sparsity structure of Jacobian:")
        P.print_nonzero()
        if(useLU==1):
            print("#  sparsity structure of LU:")
            LU.print_nonzero()

    print("#------------------------------------------------------------------#")
    print("#        OdeVec preprocessing done!                                #")
    print("#------------------------------------------------------------------#")
