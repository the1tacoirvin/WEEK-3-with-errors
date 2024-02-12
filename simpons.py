def func(x):
    # Define your function here, for example, a Gaussian-normal probability density function
    # Example: return (1 / (sigma * math.sqrt(2 * math.pi))) * math.exp(-(x - mu)**2 / (2 * sigma**2))
    pass


def simpsons_one_third_rule(func, a, b, n):
    h = (b - a) / n
    x = [a + i * h for i in range(n + 1)]
    y = [func(xi) for xi in x]

    sum_odd = sum(y[1:-1:2])
    sum_even = sum(y[2:-2:2])

    integral = (h / 3) * (y[0] + 4 * sum_odd + 2 * sum_even + y[-1])
    return integral


def main():
    a = float(input("Enter lower limit of integration: "))
    b = float(input("Enter upper limit of integration: "))
    n = int(input("Enter the number of sub-intervals (must be even): "))

    integral = simpsons_one_third_rule(func, a, b, n)
    print("Approximate integral:", integral)


if __name__ == "__main__":
    main()
