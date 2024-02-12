def luDecomposition(mat, n):
    lower = [[0] * n for _ in range(n)]
    upper = [[0] * n for _ in range(n)]

    for i in range(n):
        for k in range(i, n):
            # Compute upper triangular matrix elements
            summation = sum(lower[i][j] * upper[j][k] for j in range(i))
            upper[i][k] = mat[i][k] - summation

        for k in range(i, n):
            if i == k:
                lower[i][i] = 1  # Diagonal as 1
            else:
                # Compute lower triangular matrix elements
                summation = sum(lower[k][j] * upper[j][i] for j in range(i))
                lower[k][i] = (mat[k][i] - summation) / upper[i][i]

    # Displaying the result
    print("Lower Triangular\t\tUpper Triangular")
    for i in range(n):
        for j in range(n):
            print(lower[i][j], end="\t")
        print("\t", end="")
        for j in range(n):
            print(upper[i][j], end="\t")
        print("")


# Driver code
mat = [[1, -1, 3, 2], [-1, 5, -5, -2], [3, -5, 19, 3], [2, -2, 3, 21]]
luDecomposition(mat, 4)
