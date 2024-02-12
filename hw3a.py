def is_symmetric(matrix):
    n = len(matrix)
    for i in range(n):
        for j in range(n):
            if matrix[i][j] != matrix[j][i]:
                return False
    return True

def is_positive_definite(matrix):
    n = len(matrix)
    for i in range(n):
        if matrix[i][i] <= 0:
            return False
        for j in range(i + 1, n):
            if matrix[i][j] != matrix[j][i]:
                return False
    return True

def cholesky_decomposition(A):
    n = len(A)
    L = [[0.0] * n for _ in range(n)]

    for i in range(n):
        for j in range(i + 1):
            if i == j:
                L[i][j] = (A[i][i] - sum(L[i][k] ** 2 for k in range(j))) ** 0.5
            else:
                L[i][j] = (A[i][j] - sum(L[i][k] * L[j][k] for k in range(j))) / L[j][j]

    return L

def forward_substitution(L, b):
    n = len(L)
    y = [0.0] * n

    for i in range(n):
        y[i] = (b[i] - sum(L[i][j] * y[j] for j in range(i))) / L[i][i]

    return y

def backward_substitution(U, y):
    n = len(U)
    x = [0.0] * n

    for i in range(n - 1, -1, -1):
        x[i] = (y[i] - sum(U[i][j] * x[j] for j in range(i + 1, n))) / U[i][i]

    return x

def cholesky_solve(A, b):
    L = cholesky_decomposition(A)
    y = forward_substitution(L, b)
    x = backward_substitution([[L[j][i] for j in range(len(L))] for i in range(len(L))], y)
    return x

def doolittle_decomposition(A):
    n = len(A)
    L = [[0.0] * n for _ in range(n)]
    U = [[0.0] * n for _ in range(n)]

    for i in range(n):
        L[i][i] = 1.0
        for j in range(i, n):
            U[i][j] = A[i][j] - sum(L[i][k] * U[k][j] for k in range(i))
        for j in range(i + 1, n):
            L[j][i] = (A[j][i] - sum(L[j][k] * U[k][i] for k in range(i))) / U[i][i]

    return L, U

def solve_matrix_equation(A, b):
    if is_symmetric(A) and is_positive_definite(A):
        print("Using Cholesky Method")
        return cholesky_solve(A, b)
    else:
        print("Using Doolittle Method")
        L, U = doolittle_decomposition(A)
        y = forward_substitution(L, b)
        x = backward_substitution(U, y)
        return x

# Example matrices and vectors
A1 = [[4, 2, 2], [2, 5, 4], [2, 4, 8]]
b1 = [4, 7, 10]

A2 = [[1, 2, 3], [2, 2, 4], [3, 4, 6]]
b2 = [1, 2, 3]

# Solve the matrix equations
solution1 = solve_matrix_equation(A1, b1)
solution2 = solve_matrix_equation(A2, b2)

print("Solution for problem 1:", solution1)
print("Solution for problem 2:", solution2)
