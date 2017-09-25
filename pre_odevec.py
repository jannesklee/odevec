#!/usr/bin/env python

from sympy import symbols, Matrix, eye, fcode, SparseMatrix, lambdify

# sympy part, calculate L and U
# TODO need to include here a function which gets the rhs by krome

neq = 3
y = symbols('y(i\,1:%d)'%(neq+1))

# robertson's test
rhs = Matrix([-0.04*y[0]+1e4*y[1]*y[2],0.04*y[0]-3e7*y[1]*y[1]-1e4*y[1]*y[2],3e7*y[1]*y[1]])

jac = SparseMatrix(rhs.jacobian(y))
dt, beta = symbols('dt beta')

P = SparseMatrix(eye(neq) - dt * beta * jac)

L,U,_ = P.LUdecomposition()
# print(U.nnz())
# print(L.RL)


# Both matrices need to stored in a good proper
# CompactStorageScheme()

# write changes to file
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
    elif(srow == "#ODEVEC_RHS"):
        for i in range(rhs.shape[0]):
            fout.write("      rhs(i," + str(i + 1) + ") = " +
                       fcode(rhs[i], source_format='free', standard=95) + "\n")
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
