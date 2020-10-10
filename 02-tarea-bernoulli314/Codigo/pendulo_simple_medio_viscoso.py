"""

"""

import numpy as np 
import matplotlib.pyplot as plt 

g = 9.8 # m/s^2
l = 5.533 # m
gamma = 2.533 # s^-1


phi_0 = np.pi / gamma
omega_0 = 0

def pendulo_viscoso(t,y):
    output = np.array([y[1], -gamma * y[1] - (g/l)* np.sin(y[0])])
    return output

def k1(func, dt, t_n, y_n):
    output = dt * func(t_n, y_n)
    return output

def k2(func, dt, t_n, y_n):
    k1_n = k1(func, dt, t_n, y_n)
    k2_n = dt * func(t_n + dt/2, y_n + k1_n / 2)
    return k2_n

def k3(func, dt, t_n, y_n):
    k2_n = k2(func, dt, t_n, y_n)
    k3_n = dt * func(t_n + dt/2, y_n + k2_n / 2)
    return k3_n

def k4(func, dt, t_n, y_n):
    k3_n = k3(func, dt, t_n, y_n)
    k4_n = dt * func(t_n + dt, y_n + k3_n)
    return k4_n

def paso_rk4(func, dt, t_n, y_n):
    k1_n = k1(func, dt, t_n, y_n)
    k2_n = k2(func, dt, t_n, y_n)
    k3_n = k3(func, dt, t_n, y_n)
    k4_n = k4(func, dt, t_n, y_n)
    output = y_n + (1/6) * (k1_n + 2 * k2_n + 2 * k3_n + k4_n)
    return output

T = 5*np.pi 

dt = 0.001
t_eval_rk4 = np.arange(0, 2 * T, dt)
y_rk4 = np.zeros((len(t_eval_rk4), 2))

# cond inicial
y_rk4[0] = [phi_0, omega_0]
for i in range(1, len(t_eval_rk4)):
    y_rk4[i] = paso_rk4(pendulo_viscoso, dt, t_eval_rk4[i-1], y_rk4[i-1])

plt.figure(1)
plt.clf()

plt.plot(t_eval_rk4, y_rk4[:,0], label='Sol. RK4')

plt.xlabel('tiempo [s]')
plt.ylabel(r'$\phi(t)$', fontsize=15)
plt.legend()
plt.show()









































# Graficamos la solución de pequeñas oscilaciones
# La sol. de p.o. es phi(t) = A * e^(-1/2 * t * (factor_raiz + gamma))
#                             + B e^(1/2 * t * (factor_raiz - gamma))
# Con A,B ctes dependientes de phi_0 y omega_0

"""
T = np.pi
t_to_plot = np.linspace(0 , 2 * T, 200)
factor = gamma**2 - 4 * g / l
factor_raiz = np.sqrt(factor)
B = (omega_0 + 1/2 * factor_raiz * phi_0)/factor_raiz
A = phi_0 - B
phi_pequenas_osc = A * np.exp(-1/2 * t_to_plot * (factor_raiz + gamma)) + B * np.exp(1/2 * t_to_plot * (factor_raiz - gamma))

plt.figure(1)
plt.clf()

plt.plot(t_to_plot, phi_pequenas_osc, label='Pequ. osc.')

plt.xlabel('Tiempo[s]')
plt.ylabel(r'$\phi(t)$', fontsize=15)
plt.legend()
plt.show()
"""