import math

'''this helps do L, or lower function for cholesky decomposition'''
def forward_substitution(lower, b):
    n = len(lower)
    y = [0] * n
    for i in range(n):
        y[i] = b[i]
        for j in range(i):
            y[i] -= lower[i][j] * y[j]
        y[i] /= lower[i][i]
    return y
'''this helps do L^t, or upper function for cholesky decomposition'''
def backward_substitution(upper, y):
    n = len(upper)
    x = [0] * n
    for i in range(n - 1, -1, -1):
        x[i] = y[i]
        for j in range(i + 1, n):
            x[i] -= upper[i][j] * x[j]
        x[i] /= upper[i][i]
    return x
'''the actual cholesky decomp, the part where it multiple everything to get A, b, and x'''
def Cholesky_Decomposition(matrix, b, n):
    lower = [[0 for x in range(n)] for y in range(n)]

    # Decomposing a matrix into Lower Triangular
    for i in range(n):
        for j in range(i + 1):
            sum1 = 0
            if j == i:
                for k in range(j):
                    sum1 += pow(lower[j][k], 2)
                lower[j][j] = int(math.sqrt(matrix[j][j] - sum1))
            else:
                for k in range(j):
                    sum1 += (lower[i][k] * lower[j][k])
                if lower[j][j] > 0:
                    lower[i][j] = int((matrix[i][j] - sum1) / lower[j][j])

    # Perform forward substitution (Ly = b)
    y = forward_substitution(lower, b)
    """finds y"""

    # Perform backward substitution (L^Tx = y)
    upper = [[lower[j][i] for j in range(n)] for i in range(n)]
    x = backward_substitution(upper, y)
    """finds x value"""

    print("Solution for x:", x)
    """prints the values found"""
    print("Lower Triangular\t\tTranspose")
    for i in range(n):
        # Lower Triangular
        for j in range(n):
            print(lower[i][j], end="\t")
        print("\t", end="")

        # Transpose of Lower Triangular
        for j in range(n):
            print(lower[j][i], end="\t")
        print("")

    # Perform matrix multiplication for lower * lower^T to get the result
    result = [[sum(lower[i][k] * lower[j][k] for k in range(n)) for j in range(n)] for i in range(n)]

    print("\nResult:")
    for row in result:
        print(row)

"""prints the values found"""
print("solution to Cholesky decomposition part a:")
matrix_a = [[1, -1, 3, 2], [-1, 5, -5, -2], [3, -5, 19, 3], [2, -2, 3, 21]]
b_a = [15, -35, 94, 1]
Cholesky_Decomposition(matrix_a, b_a, 4)

print("\nsolution to Cholesky decomposition part b:")
matrix_b = [[4,2,4,0],[2,2,3,2],[4,3,6,3],[0,2,3,9]]
b_b = [20,36,60,122]
Cholesky_Decomposition(matrix_b, b_b, 4)
