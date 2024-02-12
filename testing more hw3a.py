def forward_substitution(lower, b):
    n = len(lower)
    y = [0] * n
    for i in range(n):
        y[i] = b[i]
        for j in range(i):
            y[i] -= lower[i][j] * y[j]
        y[i] /= lower[i][i]
    return y

def backward_substitution(upper, y):
    n = len(upper)
    x = [0] * n
    for i in range(n - 1, -1, -1):
        x[i] = y[i]
        for j in range(i + 1, n):
            x[i] -= upper[i][j] * x[j]
        x[i] /= upper[i][i]
    return x

def luDecomposition(mat, b):
    n = len(mat)
    lower = [[0] * n for _ in range(n)]
    upper = [[0] * n for _ in range(n)]

    for i in range(n):
        for k in range(i, n):
            summation = sum(lower[i][j] * upper[j][k] for j in range(i))
            upper[i][k] = mat[i][k] - summation

        for k in range(i, n):
            if i == k:
                lower[i][i] = 1
            else:
                summation = sum(lower[k][j] * upper[j][i] for j in range(i))
                lower[k][i] = (mat[k][i] - summation) / upper[i][i]

    # Perform forward substitution (Ly = b)
    y = forward_substitution(lower, b)

    # Perform backward substitution (Ux = y)
    x = backward_substitution(upper, y)
    print("Lower Triangular\t\tUpper Triangular")
    for i in range(n):
        for j in range(n):
            print(lower[i][j], end="\t")
        print("\t", end="")
        for j in range(n):
            print(upper[i][j], end="\t")
        print("")

    return x

# Driver code
mat = [[1, -1, 3, 2], [-1, 5, -5, -2], [3, -5, 19, 3], [2, -2, 3, 21]]
b = [15,-35,94,1]
x = luDecomposition(mat, b)
print("Solution for x:", x)