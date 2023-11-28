import numpy as np
import math
import matplotlib.pyplot as plt

#Consider the initial value problem
# y' = Ay
# y(t_0) = y_0

#y_n+1 = y_n + h(f(y_n))


#1.0
def exact_solution(y0, y, A):
    return y0 * math.exp(A * y)

#1.1
def eulerstep(A, u_old, h):
    u_new = u_old + h * A * u_old
    return u_new



def eulerint(A, y0, t0, tf, N):
    """Approximate solution using Euler method
    """
    h = (tf - t0) / N
    
    t_grid = [0] * (N+1)
    error = [0] * (N+1)
    
    t_grid[0] = y0
    
    for i in range(1, N+1):
        t_grid[i] = eulerstep(A, t_grid[i-1], h)
        error[i] = np.linalg.norm(t_grid[i] - exact_solution(y0, h * i, A))
    
    approx = t_grid[-1]
    end_err = exact_solution(y0, tf, A) - t_grid[-1]
    return t_grid, approx, end_err

print(eulerint(2,1, 0, 1, 10))
tgrid, approx, err = eulerint(2, 1, 0, 1, 10)

x_vals = np.linspace(0,1,11)
exact_vals = np.linspace(0,1,11)
for i in range(11):
    exact_vals[i] = exact_solution(1,i/10,2)
    
plt.plot(x_vals,tgrid, label ='Approximate values')
plt.plot(x_vals, exact_vals, label = 'Exact values')
plt.show()

def errVSH(A, y0, t0, tf):
    size = 100
    errors = [0] * size
    h = [0] * size
    for N in range(1,size):
        tgrid, approx, err = eulerint(A,y0,t0,tf,N)
        errors[N] = err
        h[N] = (tf - t0) / N
    print(h)
    plt.loglog(h[1:], errors[1:])
    plt.show()

errVSH(2, 1, 0, 1)