"""

"""

import numpy as np 
import matplotlib.pyplot as plt 
from scipy.integrate import solve_ivp
from mpl_toolkits.mplot3d import Axes3D

sigma = 10
beta = 8/3
rho = 28

def lorenz(s,w):
    output = np.array([sigma * (w[1] - w[0]) , w[0] * (rho - w[2]) - w[1] , w[0] * w[1] - beta * w[2]])
    return output

w_0 = np.array([1,1,1])
t_span = np.array([0.0, 10.0])
times = np.linspace(t_span[0], t_span[1], 5)

rk4 = solve_ivp(lorenz, t_span, w_0, t_eval=times)

fig = plt.figure(1)
fig.clf()

ax = fig.add_subplot(111, projection='3d')

ax.plot(rk4.y[0], rk4.y[1], rk4.y[2])

w_02 = np.array([1,1.0533,1])

rk4_2 = solve_ivp(lorenz, t_span, w_02, t_eval=times)


ax.plot(rk4_2.y[0], rk4_2.y[1], rk4_2.y[2])

ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')

plt.show()