"""
Este script integra la ecuacion de movimiento de un pendulo real
usando el algoritmo de Runge-Kutta de orden 2. 
"""

import numpy as np
import matplotlib.pyplot as plt

g = 9.8  # en m/s^2, aceleracion de gravedad
R = 1. # en m, largo del pendulo

phi_0 = np.pi / 2
omega_0 = 10

# Graficamos la solucion de pequeñas oscilaciones

T = 2 * np.pi / np.sqrt(g / R)  # periodo pequeñas oscilaciones
t_to_plot = np.linspace(0, 2 * T, 200)
phi_pequenas_osc = phi_0 * np.cos(np.sqrt(g/R) * t_to_plot)


# Implementando solucion usando RK2

# dy/dt = f(t, y)
def func_pendulo(t, y):
    """ ec. del pendulo real
    d2t/dt2 = -g/R sin(y[0])

    IMPUTS:
    =======
    t : [float], tiempo
    y : [np.ndarray], [phi, dphi/dt]
    """
    output = np.array([y[1], -g / R * np.sin(y[0])])
    return output

def K1(func, dt, t_n, y_n):
    """constante K1 del algoritmo de runge kutta 2
    """
    output = dt * func(t_n, y_n)
    return output


def K2(func, dt, t_n, y_n):
    """ constante K2 del algoritmo de runge kutta 2
    """
    k1_n = K1(func,dt, t_n, y_n)
    k2_n = dt * func(t_n + dt/2, y_n + k1_n / 2)
    return k2_n


def paso_rk2(func, dt, t_n, y_n):
    output = y_n + K2(func, dt, t_n, y_n)
    return output


dt = 0.001
t_eval_rk2 = np.arange(0, 2 * T, dt)
y_rk2 = np.zeros((len(t_eval_rk2), 2))

# cond inicial
y_rk2[0] = [phi_0, omega_0]
for i in range(1, len(t_eval_rk2)):
    y_rk2[i] = paso_rk2(func_pendulo, dt, t_eval_rk2[i-1], y_rk2[i-1])


plt.figure(1)
plt.clf()

# plt.plot(t_to_plot, phi_pequenas_osc, label='peq osc')
plt.plot(t_eval_rk2, y_rk2[:,0], label='sol rk2')

plt.xlabel('tiempo [s]')
plt.ylabel(r'$\phi(t)$', fontsize=15)
plt.legend()
plt.show()
plt.savefig('pequenas_osc_vs_pendulo_real.png')