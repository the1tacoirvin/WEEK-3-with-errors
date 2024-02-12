def doolittle_lu_decomposition(A):
    n = len(A)
    L = [[0.0] * n for _ in range(n)]
    U = [[0.0] * n for _ in range(n)]

    for i in range(n):
        # Upper triangular matrix
        for k in range(i, n):
            sum_val = 0
            for j in range(i):
                sum_val += L[i][j] * U[j][k]
            U[i][k] = A[i][k] - sum_val

        # Lower triangular matrix
        for k in range(i, n):
            if i == k:
                L[i][i] = 1.0
            else:
                sum_val = 0
                for j in range(i):
                    sum_val += L[k][j] * U[j][i]
                L[k][i] = (A[k][i] - sum_val) / U[i][i]

    return L, U
def doolittle_solve(A, b):
    L, U = doolittle_lu_decomposition(A)

    # Solve Ly = b using forward substitution
    n = len(b)
    y = [0.0] * n
    for i in range(n):
        sum_val = 0
        for j in range(i):
            sum_val += L[i][j] * y[j]
        y[i] = b[i] - sum_val

    # Solve Ux = y using backward substitution
    x = [0.0] * n
    for i in range(n-1, -1, -1):
        sum_val = 0
        for j in range(i+1, n):
            sum_val += U[i][j] * x[j]
        x[i] = (y[i] - sum_val) / U[i][i]

    return x

