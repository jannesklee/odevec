import sympy as sym

def load_rober():
    print("#  Setup: ROBER test")
    neq = 3
    k = 0
    nrea = 0

    # define sympy symbols
    y = list(sym.symbols('y(\:\,1:%d)'%(neq+1)))

    rhs = sym.Matrix(
            [-0.04*y[0]+1e4*y[1]*y[2],
             0.04*y[0]-3e7*y[1]*y[1]-1e4*y[1]*y[2],
             3e7*y[1]*y[1]])

    return (neq, k, nrea, y, rhs)


def load_orego():
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

    return (neq, k, nrea, y, rhs)


def load_primordial():
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

    return (neq, k, nrea, y, rhs)
