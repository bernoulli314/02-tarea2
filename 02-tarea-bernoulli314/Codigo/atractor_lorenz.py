"""
PROBLEMA 2
Algoritmo que soluciona un caso del conocido atractor de Lorenz
cuyas ecuaciones son
dx/ds = sigma(y-x)
dy/ds = x(rho-z)
dz/ds = xy - beta*z
En este caso se considerarán los parámetros para obtener la solución mas famosa
sigma = 10
beta = 8/3
rho = 28
Esto se resolverá mediante el método RK4, pero esta vez desde la
libreria scipy.integrate
"""
# Diego Gonzalez Flores
# RUT: 20.300.533-4

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from mpl_toolkits.mplot3d import Axes3D

# se definen las constantes
sigma = 10
beta = 8/3
rho = 28


def lorenz(s, w):
    """
    Ecuación de atractor de Lorenz mediante arreglos de numpy con
    formato necesario para implementar RK4
    PARAMETROS
    s: tiempo
    w: arreglo numpy del tipo np.array([x,y,z])
    """
    output = np.array([sigma * (w[1] - w[0]),
                      w[0] * (rho - w[2]) - w[1], w[0] * w[1] - beta * w[2]])
    return output

w_0 = np.array([1, 1, 1])           # condiciones iniciales
s_span = np.array([0.0, 100.0])         # valor inicial y final de s
times = np.linspace(s_span[0], s_span[1], 5000)        # cant de evaluaciones

rk4 = solve_ivp(lorenz, s_span, w_0, t_eval=times)      # rk4

fig = plt.figure(1)
fig.clf()

ax = fig.add_subplot(111, projection='3d')

ax.plot(rk4.y[0], rk4.y[1], rk4.y[2])       # plot c.i.

w_02 = np.array([1, 1.0533, 1])         # condiciones iniciales perturbado

rk4_2 = solve_ivp(lorenz, s_span, w_02, t_eval=times)       # rk4


ax.plot(rk4_2.y[0], rk4_2.y[1], rk4_2.y[2])         # plot c.i. perturbada

ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')

plt.show()

plt.figure(2)
plt.clf()

plt.plot(rk4.y[1], rk4.y[2], label='Condiciones Iniciales')
plt.plot(rk4_2.y[1], rk4_2.y[2], label='Condiciones Perturbadas')
plt.xlabel('Coordenada Y')
plt.ylabel('Coordenada Z')
plt.legend()
plt.show()