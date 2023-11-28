
import numpy as np
import matplotlib.pyplot as plt
import math

from cProfile import label
import numpy as np
import matplotlib.pyplot as plt
import math

a = 2

def f(t, y):
    #y'(t) = ay
    return a * y

#y(t) = y_0*e^(at)
def actual(t_n, y0):
    return y0 * math.exp(a*t_n)
    

#f(t0) = y0
#COMPUTE FOUR DIFFERENT DERIVATIVES to form average derivative

def newstep(tol, curr_r, prev_r, prev_h, k):
    return (tol / curr_r)**(2 / (3 * k)) * (tol / prev_r)**(-1/(3*k)) * prev_h

def RK34Step(func, t_n, y_n, h):
    
    Y1 = func(t_n, y_n)
    Y2 = func(t_n + h / 2, y_n + h * Y1 / 2)
    Y3 = func(t_n + h / 2, y_n + h * Y2 / 2)
    Y4 = func(t_n + h, y_n + h* Y3)

    Z3 = func(t_n + h, y_n - h*Y1 + 2*h*Y2)
    
    y_next = y_n + (h/6) * (Y1 + 2 * Y2 + 2 * Y3 + Y4)
    l_next = h / 6 * (2*Y2 + Z3 - 2*Y3 - Y4)
    
    return [y_next, l_next]

def adaptiveRK34(f, t0, tf, y0, tol):
    h = abs(tf - t0)*tol**(1/4) / (100 * (1 + abs(f(t0, y0))))   


    T = [t0]
    Y = [y0]

    prev_r = tol
    
    k = 4   

    while T[-1] < tf:
        t = min(T[-1] + h, tf)
        h = t - T[-1]
        y_n, l_n = RK34Step(f, t, Y[-1], h)
        r_n = max(abs(l_n), 0.000001) #Floating point error
    
        T.append(t)
        Y.append(y_n)
        h = newstep(tol, r_n, prev_r, h, k)
        prev_r = r_n

    return [T, Y]

tol = 100
t0 = 0
tf = 1
y0 = 1

t, y = adaptiveRK34(f, t0, tf, y0, tol)

actual_t = np.linspace(t0, tf, 1000)
actual_y = [actual(t, y0) for t in actual_t]

plt.plot(t, y, label="approx y")
plt.plot(actual_t, actual_y, label="actual y")
plt.title('TOL = 100')
plt.ylabel('y')
plt.xlabel('x')
plt.legend()
plt.show()


#Task 1.1

def part11():
    
    def analytical_solution(t, y0, lambda_val):
        return y0 * math.exp(lambda_val * t)

    a = 2
    def f(t, y):
        #y'(t) = ay
        return a * y
    
    def RK4step(func, t_n, y_n, h):
    
        Y1 = func(t_n, y_n)
        Y2 = func(t_n + h / 2, y_n + h * Y1 / 2)
        Y3 = func(t_n + h / 2, y_n + h * Y2 / 2)
        Y4 = func(t_n + h, y_n + h* Y3)
        
        y_next = y_n + (h/6) * (Y1 + 2 * Y2 + 2 * Y3 + Y4)
        
        err_est = (h/6) * (2 * Y2)
        
        return y_next

    def err_on_h():
        t0, t_final = 0, 1
        y0 = 1
        lambda_val = 2
        n_values = np.logspace(1,3,100)
        n_values = [int(n) for n in n_values]
        h_values = [1/n for n in n_values]
        err = []

        for n in n_values:
            y_num = y0
            t = 0
            h = 1 / n
            for _ in range(n):
                y_num = RK4step(f, t, y_num, h)
                t += h
                
            y_analytical = analytical_solution(t_final, y0, lambda_val)
            error = abs(y_num - y_analytical)
            err.append(error)
            
        h4 = [a ** 4 for a in h_values]
        h3 = [a ** 3 for a in h_values]
        h5 = [a ** 5 for a in h_values]
            
        plt.loglog(h_values, err, marker='o', label='Global Error')
        plt.loglog(h_values, h3, marker='o', label='h3')
        plt.loglog(h_values, h4, marker='o', label='h4')
        plt.loglog(h_values, h5, marker='o', label='h5')
        plt.xlabel('Step size (h)')
        plt.ylabel('Error')
        plt.title('Global Error of RK4 Method')
        plt.legend()
        plt.show()
    err_on_h()
    
part11()

a = 2
y_0 = 1

def f(t, y):
    #y'(t) = ay
    return a * y

#y(t) = y_0*e^(at)
def actual(t_n):
    return y_0 * math.exp(a*t_n)
    

#f(t0) = y0
#COMPUTE FOUR DIFFERENT DERIVATIVES to form average derivative

def newstep(tol, curr_r, prev_r, prev_h, k):
    return (tol / curr_r)**(2 / (3 * k)) * (tol/prev_r)**(-1/(3*k)) * prev_hold

def RK34step(func, t_n, y_n, h):
    
    Y1 = func(t_n, y_n)
    Y2 = func(t_n + h / 2, y_n + h * Y1 / 2)
    Y3 = func(t_n + h / 2, y_n + h * Y2 / 2)
    Y4 = func(t_n + h, y_n + h* Y3)

    Z3 = func(t_n + h, y_n - h*Y1 + 2*h*Y2)
    
    y_next = y_n + (h/6) * (Y1 + 2 * Y2 + 2 * Y3 + Y4)
    l_next = h / 6 * (2*Y2+Z3 - 2*Y3 - Y4)
    
    return [y_next, l_next]

def adaptiveRK34(f, t0, tf, y0, tol):
    h = abs(tf - t0)*tol**(1/4) / (100 * (1 + abs(f(y_0))))   
    T = [t0]
    Y = [y_0]

    prev_r = 0
    
    k = 4

    while T[-1] < tf:
        t = min(T[-1] + h, tf)
        h = t - T[-1]
        y_n, l_n = RK34Step(f, t, Y[-1], h)
        r_n = abs(l_n)
    
        T.append(t)
        Y.append(y_n)
         
        h = newstep(tol, r_n, prev_r, h, k)
        prev_r = r_n
    
    return [T, Y]
    
         



    

def endpoint_err(start, stop, steps):
    T = np.linspace(start,stop,steps)
    h = T[1] - T[0]
    approx_y = [y_0]
    actual_y = [y_0]
    for t in T[1:]:
        approx_y.append(RK4step(f, t, approx_y[-1], h))
        actual_y.append(actual(t))

    diff = [abs(approx_y[i] - actual_y[i]) for i in range(len(T))]
    return diff[-1], h #endpoint error


            
            
            
            
            







