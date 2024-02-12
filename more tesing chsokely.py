import math


def Cholesky_Decomposition(matrix, n):
    lower = [[0 for x in range(n)]
             for y in range(n)]

    # Decomposing a matrix
    # into Lower Triangular
    for i in range(n):
        for j in range(i + 1):
            sum1 = 0

            # summation for diagonals
            if (j == i):
                for k in range(j):
                    sum1 += pow(lower[j][k], 2)
                lower[j][j] = int(math.sqrt(matrix[j][j] - sum1))
            else:
                # Evaluating L(i, j)
                # using L(j, j)
                for k in range(j):
                    sum1 += (lower[i][k] * lower[j][k])
                if (lower[j][j] > 0):
                    lower[i][j] = int((matrix[i][j] - sum1) / lower[j][j])

    # Displaying Lower Triangular
    # and its Transpose
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


# Driver Code
n = 4
matrix = [[1, -1, 3, 2], [-1, 5, -5, -2], [3, -5, 19, 3], [2, -2, 3, 21]]
b = [15,-35,94,1]
Cholesky_Decomposition(matrix, n)
