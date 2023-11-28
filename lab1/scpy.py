from scipy.integrate import solve_ivp as sol
import matplotlib.pyplot as plt

mys = [ 10, 15, 22, 33, 47, 68, 100, 150, 220, 330, 470, 680, 1000]

times = []
for i in range(len(mys)):
    my = mys[i]
    t0 = 0
    tf = 0.7 * my
    y0 = [2, 0]

    def dydt(t, y):
        x, y = y
        l = y
        r = my * (1-x**2) * y - x
        return [l, r]

    ans = sol(dydt, (t0, tf), y0, method='BDF')
    #print(ans)
    #print(len(ans.t))
    times.append(len(ans.t))
    
plt.plot(mys, times)

plt.title('time and my')
plt.ylabel('time')
plt.xlabel('my')
plt.show()

    