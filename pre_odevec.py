#!/usr/bin/env python

from sympy import symbols, Matrix, eye, fcode, SparseMatrix, lambdify, nsimplify
import numpy as np

# sympy part, calculate L and U
# TODO need to include here a function which gets the rhs by krome


#------------------------------------------------------------------------------#
# robertson's test
#neq = 3
#y = symbols('y(i\,1:%d)'%(neq+2))
#
#rhs = Matrix([-0.04*y[0]+1e4*y[1]*y[2],0.04*y[0]-3e7*y[1]*y[1]-1e4*y[1]*y[2],3e7*y[1]*y[1]])
#------------------------------------------------------------------------------#
# primordial cooling

nvector = 1
neq = 12
nrea = 38
maxorder = 5

idx_E=1
idx_Hk=2
idx_H=3
idx_HE=4
idx_H2=5
idx_D=6
idx_HD=7
idx_Hj=8
idx_HEj=9
idx_H2j=10
idx_Dj=11
idx_HEjj=12
idx_CR=13
idx_g=14
idx_Tgas=15
idx_dummy=16

y = symbols('y(i\,0:%d)'%(neq+1))
k = symbols('k(i\,0:%d)'%(nrea+1))

rhs = Matrix([
    -k[1]*y[idx_H]*y[idx_E] +2.e0*k[1]*y[idx_H]*y[idx_E] -k[2]*y[idx_Hj]*y[idx_E]-k[3]*y[idx_Hj]*y[idx_E] -k[4]*y[idx_HE]*y[idx_E] +2.e0*k[4]*y[idx_HE]*y[idx_E] -k[5]*y[idx_HEj]*y[idx_E] -k[6]*y[idx_HEj]*y[idx_E] -k[7]*y[idx_HEj]*y[idx_E] +2.e0*k[7]*y[idx_HEj]*y[idx_E] -k[8]*y[idx_HEjj]*y[idx_E] -k[9]*y[idx_H]*y[idx_E] +k[10]*y[idx_Hk]*y[idx_H] +k[11]*y[idx_Hk]*y[idx_H] -k[16]*y[idx_H2]*y[idx_E] +k[16]*y[idx_H2]*y[idx_E] -k[18]*y[idx_Hk]*y[idx_E] +2.e0*k[18]*y[idx_Hk]*y[idx_E] +k[19]*y[idx_Hk]*y[idx_H] +k[20]*y[idx_Hk]*y[idx_H] +k[22]*y[idx_Hk]*y[idx_Hj] -k[23]*y[idx_H2j]*y[idx_E] -k[24]*y[idx_H2j]*y[idx_E] +k[37]*y[idx_D]*y[idx_Hk] -k[38]*y[idx_Dj]*y[idx_E],
    +k[9]*y[idx_H]*y[idx_E] -k[10]*y[idx_Hk]*y[idx_H] -k[11]*y[idx_Hk]*y[idx_H] -k[18]*y[idx_Hk]*y[idx_E] -k[19]*y[idx_Hk]*y[idx_H] -k[20]*y[idx_Hk]*y[idx_H] -k[21]*y[idx_Hk]*y[idx_Hj] -k[22]*y[idx_Hk]*y[idx_Hj] -k[25]*y[idx_H2j]*y[idx_Hk] -k[37]*y[idx_D]*y[idx_Hk],
    -k[1]*y[idx_H]*y[idx_E] +k[2]*y[idx_Hj]*y[idx_E] +k[3]*y[idx_Hj]*y[idx_E] -k[9]*y[idx_H]*y[idx_E] -k[10]*y[idx_Hk]*y[idx_H] -k[11]*y[idx_Hk]*y[idx_H] -k[12]*y[idx_H]*y[idx_Hj] -k[13]*y[idx_H]*y[idx_Hj] -k[14]*y[idx_H2j]*y[idx_H] +k[15]*y[idx_H2]*y[idx_Hj] +2.e0*k[16]*y[idx_H2]*y[idx_E] -k[17]*y[idx_H2]*y[idx_H] +3.e0*k[17]*y[idx_H2]*y[idx_H] +k[18]*y[idx_Hk]*y[idx_E] -k[19]*y[idx_Hk]*y[idx_H] +2.e0*k[19]*y[idx_Hk]*y[idx_H] -k[20]*y[idx_Hk]*y[idx_H] +2.e0*k[20]*y[idx_Hk]*y[idx_H] +2.e0*k[21]*y[idx_Hk]*y[idx_Hj] +2.e0*k[23]*y[idx_H2j]*y[idx_E] +2.e0*k[24]*y[idx_H2j]*y[idx_E] +k[25]*y[idx_H2j]*y[idx_Hk] -3.e0*k[26]*y[idx_H]*y[idx_H]*y[idx_H] +k[26]*y[idx_H]*y[idx_H]*y[idx_H] -3.e0*k[27]*y[idx_H]*y[idx_H]*y[idx_H] +k[27]*y[idx_H]*y[idx_H]*y[idx_H] -2.e0*k[28]*y[idx_H2]*y[idx_H]*y[idx_H] -2.e0*k[29]*y[idx_H2]*y[idx_H]*y[idx_H] +k[30]*y[idx_Hj]*y[idx_D] -k[31]*y[idx_H]*y[idx_Dj] +k[34]*y[idx_H2]*y[idx_D] +k[35]*y[idx_H2]*y[idx_D] -k[36]*y[idx_HD]*y[idx_H],
    -k[4]*y[idx_HE]*y[idx_E] +k[5]*y[idx_HEj]*y[idx_E] +k[6]*y[idx_HEj]*y[idx_E],
    +k[10]*y[idx_Hk]*y[idx_H] +k[11]*y[idx_Hk]*y[idx_H] +k[14]*y[idx_H2j]*y[idx_H] -k[15]*y[idx_H2]*y[idx_Hj] -k[16]*y[idx_H2]*y[idx_E] -k[17]*y[idx_H2]*y[idx_H] +k[25]*y[idx_H2j]*y[idx_Hk] +k[26]*y[idx_H]*y[idx_H]*y[idx_H] +k[27]*y[idx_H]*y[idx_H]*y[idx_H] -k[28]*y[idx_H2]*y[idx_H]*y[idx_H] +2.e0*k[28]*y[idx_H2]*y[idx_H]*y[idx_H] -k[29]*y[idx_H2]*y[idx_H]*y[idx_H] +2.e0*k[29]*y[idx_H2]*y[idx_H]*y[idx_H] -k[32]*y[idx_H2]*y[idx_Dj] +k[33]*y[idx_HD]*y[idx_Hj] -k[34]*y[idx_H2]*y[idx_D] -k[35]*y[idx_H2]*y[idx_D] +k[36]*y[idx_HD]*y[idx_H],
    -k[30]*y[idx_Hj]*y[idx_D] +k[31]*y[idx_H]*y[idx_Dj] -k[34]*y[idx_H2]*y[idx_D] -k[35]*y[idx_H2]*y[idx_D] +k[36]*y[idx_HD]*y[idx_H] -k[37]*y[idx_D]*y[idx_Hk] +k[38]*y[idx_Dj]*y[idx_E],
    +k[32]*y[idx_H2]*y[idx_Dj] -k[33]*y[idx_HD]*y[idx_Hj] +k[34]*y[idx_H2]*y[idx_D] +k[35]*y[idx_H2]*y[idx_D] -k[36]*y[idx_HD]*y[idx_H] +k[37]*y[idx_D]*y[idx_Hk],
    +k[1]*y[idx_H]*y[idx_E] -k[2]*y[idx_Hj]*y[idx_E] -k[3]*y[idx_Hj]*y[idx_E] -k[12]*y[idx_H]*y[idx_Hj] -k[13]*y[idx_H]*y[idx_Hj] +k[14]*y[idx_H2j]*y[idx_H] -k[15]*y[idx_H2]*y[idx_Hj] -k[21]*y[idx_Hk]*y[idx_Hj] -k[22]*y[idx_Hk]*y[idx_Hj] -k[30]*y[idx_Hj]*y[idx_D] +k[31]*y[idx_H]*y[idx_Dj] +k[32]*y[idx_H2]*y[idx_Dj] -k[33]*y[idx_HD]*y[idx_Hj],
    +k[4]*y[idx_HE]*y[idx_E] -k[5]*y[idx_HEj]*y[idx_E] -k[6]*y[idx_HEj]*y[idx_E] -k[7]*y[idx_HEj]*y[idx_E] +k[8]*y[idx_HEjj]*y[idx_E],
    +k[12]*y[idx_H]*y[idx_Hj] +k[13]*y[idx_H]*y[idx_Hj] -k[14]*y[idx_H2j]*y[idx_H] +k[15]*y[idx_H2]*y[idx_Hj] +k[22]*y[idx_Hk]*y[idx_Hj] -k[23]*y[idx_H2j]*y[idx_E] -k[24]*y[idx_H2j]*y[idx_E] -k[25]*y[idx_H2j]*y[idx_Hk],
    +k[30]*y[idx_Hj]*y[idx_D] -k[31]*y[idx_H]*y[idx_Dj] -k[32]*y[idx_H2]*y[idx_Dj] +k[33]*y[idx_HD]*y[idx_Hj] -k[38]*y[idx_Dj]*y[idx_E],
    +k[7]*y[idx_HEj]*y[idx_E] -k[8]*y[idx_HEjj]*y[idx_E]
    ])


#------------------------------------------------------------------------------#

jac = nsimplify(SparseMatrix(rhs.jacobian(y)[:,1:]))
dt, beta = symbols('dt beta')

P = SparseMatrix(eye(neq) - dt * beta * jac)


if(np.shape(P)[0] < 5):
    #L,U,_ = P.LUdecomposition()
    Piv,L,D,U = P.LUdecompositionFF()
    # print(U.nnz())
    # print(L.RL)


# Both matrices need to stored in a good proper
# CompactStorageScheme()

# write changes to file
fh_list = [open("odevec.F90"),open("odevec_commons.F90")]
fout    = [open("odevec.f90","w"),open("odevec_commons.f90","w")]

for k,fh in enumerate(fh_list):
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
            fout[k].write("      integer :: nvector=" + str(nvector) + "\n")
        elif(srow == "#ODEVEC_EQUATIONS"):
            fout[k].write("      integer :: neq=" + str(neq) + "\n")
        elif(srow == "#ODEVEC_MAXORDER"):
            fout[k].write("      integer :: maxorder=" + str(maxorder) + "\n")
        elif(srow == "#ODEVEC_REACTIONS"):
            fout[k].write("      integer :: nrea=" + str(nrea) + "\n")
        else:
            srow = row.strip()
            fout[k].write(row)


print("sparsity structure of P:")
P.print_nonzero()

if(np.shape(P)[0] < 5):
    print("sparsity structure of L:")
    L.print_nonzero()

if(np.shape(P)[0] < 5):
    print("sparsity structure of U:")
    U.print_nonzero()

print("OdeVec preprocessing done!")
