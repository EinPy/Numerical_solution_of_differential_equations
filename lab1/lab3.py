import math
import matplotlib.pyplot as plt
import sys
import time
import numpy as np
from collections import defaultdict

a, b, c, d = 3, 9, 15, 15

my = [ 10, 15, 22, 33, 47, 68, 100, 150, 220, 330, 470, 680, 1000]


def dx(x, y):
    return a * x - b * x * y

def dy(x, y):
    return c * x * y - d * y

def lotka(t, u, my):
    x, y = u
    l = y
    r = my * (1-x**2) * y -
    return np.array([l, r])

def H(x, y):
    return c * x + b * y - d * math.log(x) - a * math.log(y)

def vec_len(v):
    a, b = v
    return math.sqrt(a * a + b * b)

def stepL(f, t, u):
    return math.sqrt(f(t, u)[0] ** 2 + f(t, u)[1] ** 2)


#f(t0) = y0
#COMPUTE FOUR DIFFERENT DERIVATIVES to form average derivative
def newstep(tol, curr_r, prev_r, prev_h, k):
    return (tol / curr_r)**(2 / (3 * k)) * (tol / prev_r)**(-1/(3*k)) * prev_h

def RK34Step(func, t_n, y_n, h):
    y_n = np.array(y_n)
    Y1 = func(t_n, y_n)
    #print(Y1)
    Y2 = func(t_n + h / 2, y_n + h * Y1 / 2)
    Y3 = func(t_n + h / 2, y_n + h * Y2 / 2)
    Y4 = func(t_n + h, y_n + h* Y3)

    Z3 = func(t_n + h, y_n - h*Y1 + 2*h*Y2)
    
    y_next = y_n + (h/6) * (Y1 + 2 * Y2 + 2 * Y3 + Y4)
    l_next = h / 6 * (2*Y2 + Z3 - 2*Y3 - Y4)
    
    return [y_next, vec_len(l_next)]


def adaptiveRK34(f, t0, tf, u0, tol):
    h = abs(tf - t0)*tol**(1/4) / (100 * (1 + stepL(f, t0, u0)))   

    T = [t0]
    U = [u0]
    x, y = [u0[0]],[u0[1]]

    prev_r = tol
    
    k = 4   

    while T[-1] < tf:
        t = min(T[-1] + h, tf)
        h = t - T[-1]
        u_n, l_n = RK34Step(f, t, U[-1], h)
        #print(l_n)
        r_n = max(abs(l_n), 0.000000000001) #Floating point error
    
        T.append(t)
        U.append(u_n)
        x.append(u_n[0])
        y.append(u_n[1])
        h = newstep(tol, r_n, prev_r, h, k)
        prev_r = r_n

    return [T, x, y]

tol = 10**(-8)
t0 = 0
tf = 10


u0 = np.array([1, 1/3 + 0.1])

t, x, y = adaptiveRK34(lotka, t0, tf, u0, tol)

period = []
seen = defaultdict(None)
for i, (a, b) in enumerate(zip(x, y)):
    if abs(u0[0] - a) < 0.001 and abs(u0[1] - b) < 0.001:
        period.append((i, a, b))
        if len(period) == 4:
            break
        


print(period)
#length first hit at 258 steps
print(t[period[-1][0]])
#time 1.0226

plt.plot(t, x, label="rabbits")
plt.plot(t, y, label="foxes")
plt.title('Foxes and rabbits')
plt.ylabel('Fox / Rabbit')
plt.xlabel('time')
plt.legend()
plt.show()


plt.plot(x, y)
for _, a, b in period:
    plt.plot(a, b, color='red', marker='o')
    
plt.title('Foxes and rabbits')
plt.ylabel('foxes')
plt.xlabel('rabbits')
plt.show()





# def drift(x, y):
#     return abs(H(x, y) / H(1, 1) - 1)
    
    
# tol = 10**(-8)
# t0 = 0
# tf = 1000


# u0 = np.array([1, 1])

# t, x, y = adaptiveRK34(lotka, t0, tf, u0, tol)
# Harr = [drift(a,b) for a,b in zip(x, y)]
# plt.plot(t, Harr)
# plt.title('Diff')
# plt.ylabel('H/H_0')
# plt.xlabel('time')
# plt.show()




