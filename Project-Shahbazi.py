# Amir Shahbazi - 9812033 - Final Project

import numpy as np
import math

# ----------------------------------------------------------------------------------------
delta_t = 10**(-4)
delta_x = [0.2, 0.1, 0.05, 0.025]

alpha = 1

# ---------------------------------------------------------------------------------------
# طرح تفاضلات متناهی پيشرو

for h in delta_x:
    lambdaa = ((alpha**2)*delta_t)/(h**2)
    n = int(1/h)
    m = int(0.1/delta_t)
    A = np.zeros((n-1, n-1))

    for i in range(n-1):
        A[i][i] = 1-(2*lambdaa)
        if 0 <= i-1 <= n-2:
            A[i-1][i] = lambdaa
        if 0 <= i+1 <= n-2:
            A[i+1][i] = lambdaa

    answer = np.zeros((m+1, n+1))

    for i in range(m+1):
        answer[i][0] = 0
        answer[i][n] = 0

    x = 0
    for i in range(1, n):
        x += h
        y = (3*math.sin((math.pi)*x))+5*(math.sin((4*math.pi)*x))

        answer[0][i] = y

    w0 = np.transpose(answer)[1:n, 0]
    for i in range(1, m):
        w = np.matmul(A, w0)

        answer[i][1:n] = np.transpose(w)
        w0 = w

    print("روش تفاضلات متناهی پیشرو")
    print("---------------------------")
    print("\n")
    for i in range(m+1):
        for j in range(n+1):
            x = j*h
            t = i*delta_t

            u_x_t = 3*(np.exp(-(math.pi**2)*t))*(math.sin(math.pi*x)) + \
                5*(np.exp(-16*(math.pi**2)*t))*(math.sin(4*math.pi*x))

            wi = answer[i][j]

            error = abs(wi - u_x_t)

            print(f"x={x} , t={t} , Error={error}")

    print("\n")
    print("\n")


# ---------------------------------------------------------------------------------------
# طرح تفاضلات متناهی پسرو

for h in delta_x:
    lambdaa = ((alpha**2)*delta_t)/(h**2)
    n = int(1/h)
    m = int(0.1/delta_t)
    A = np.zeros((n-1, n-1))

    for i in range(n-1):
        A[i][i] = 1+(2*lambdaa)
        if 0 <= i-1 <= n-2:
            A[i-1][i] = -lambdaa
        if 0 <= i+1 <= n-2:
            A[i+1][i] = -lambdaa

    answer = np.zeros((m+1, n+1))

    for i in range(m+1):
        answer[i][0] = 0
        answer[i][n] = 0

    x = 0
    for i in range(1, n):
        x += h
        y = (3*math.sin((math.pi)*x))+5*(math.sin((4*math.pi)*x))

        answer[0][i] = y

    w0 = answer[0, 1:n]
    for i in range(1, m):
        w = np.linalg.solve(A, w0)

        answer[i][1:n] = w
        w0 = w

    print("---------------------------")
    print("روش تفاضلات متناهی پسرو")
    print("---------------------------")
    print("\n")

    for i in range(m+1):
        for j in range(n+1):
            x = j*h
            t = i*delta_t

            u_x_t = 3*(np.exp(-(math.pi**2)*t))*(math.sin(math.pi*x)) + \
                5*(np.exp(-16*(math.pi**2)*t))*(math.sin(4*math.pi*x))

            wi = answer[i][j]

            error = abs(wi - u_x_t)

            print(f"x={x} , t={t} , Error={error}")

    print("\n")
    print("\n")
