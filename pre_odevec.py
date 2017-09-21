#!/usr/bin/env python

from shutil import copyfile
from sympy import symbols, Matrix, eye, fcode, SparseMatrix, lambdify

# sympy part, calculate L and U
# TODO need to include here a function which gets the jacobian by
# KROME

neq = 3
x = symbols('y(i\,0:%d)' % (neq + 1))
dt, beta = symbols('dt beta')
J = SparseMatrix([[-4 / 100, 10**4 * x[3], 10**4 * x[2]],
                  [4 / 100, -6 * 10**7 * x[2] - 10**4 * x[3], -10**4 * x[2]],
                  [0, 6 * 10**7 * x[2], 0]])

P = SparseMatrix(eye(neq) - dt * beta * J)

L,U,_ = P.LUdecomposition()
# print(U.nnz())
# print(L.RL)


# Both matrices need to stored in a good proper
# CompactStorageScheme()

# write changes to file
copyfile("bdf_NEC.F90", "bdf_NEC.f90")

fh = open("bdf_NEC.F90")
fout = open("bdf_NEC.f90", "w")

for row in fh:
    srow = row.strip()
    if(srow == "#ODEVEC_L"):
        for i in range(L.shape[0]):
            for j in range(L.shape[1]):
                fout.write("      L(i," + str(i + 1) + "," + str(j + 1) + ") = " +
                           fcode(L[i, j], source_format='free', standard=95) + "\n")
    elif(srow == "#ODEVEC_U"):
        for i in range(L.shape[0]):
            for j in range(L.shape[1]):
                fout.write("      U(i," + str(i + 1) + "," + str(j + 1) + ") = " +
                           fcode(U[i, j], source_format='free', standard=95) + "\n")
    else:
        srow = row.strip()
        fout.write(row)


print("sparsity structure of P:")
P.print_nonzero()

print("sparsity structure of L:")
L.print_nonzero()

print("sparsity structure of U:")
U.print_nonzero()

print("OdeVec preprocessing done!")
