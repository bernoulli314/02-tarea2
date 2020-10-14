"""
PROBLEMA 1
Algoritmo que soluciona el problema de un péndulo simple
en un medio viscoso con ecuación y'' = -ay' - (g / l) sin(y)
(a partir de ahora phi) mediante el método RK4,
el cual es programado en el mismo script.
"""
# Diego Gonzalez Flores
# RUT: 20.300.533-4

import numpy as np
import matplotlib.pyplot as plt


g = 9.8         # m/s^2
l = 5.533       # m
gamma = 2.533       # s^-1


phi_0 = np.pi / gamma       # phi inicial
phipunto_0 = 0          # dphi/dt inicial

# Graficamos la solución de pequeñas oscilaciones
# Su ecuación es phi'' = - gamma * phi' - (g/l) phi
# La sol. de p.o. es phi(t) = A * e^(-1/2 * t * gamma) * cos(omega * t + phi_0)
# Con A cte dependiente de phi_0

omega = np.sqrt(g / l - (gamma ** 2) / 4)
A = phi_0 / np.cos(phi_0)
T = 2 * np.pi / omega       # Periodo

t_to_plot = np.linspace(0, 2 * T, 200)

phi_pequenas_osc = A * np.e ** (- t_to_plot * gamma / 2) \
                  * np.cos(omega * t_to_plot + phi_0)


def pendulo_viscoso(t, y):
    """
    Ecuación péndulo en medio viscoso mediante arreglos de numpy con
    formato necesario para implementar RK4
    PARAMETROS
    t: tiempo
    y: arreglo numpy del tipo np.array([a,b]) con b = da/dt
    """
    output = np.array([y[1], -gamma * y[1] - (g / l) * np.sin(y[0])])
    return output


# Se definen los pasos necesarios para implementar el método RK4
# Cuya ecuación es
# y_n+1 = y_i + (1/6) * h * (k_1 + 2*k_2 + 2_k_3 + k_4)

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

T = 2*np.pi / omega         # Periodo(se usa el mismo creado anteriormente)

dt = 0.0001         # dt = h
t_eval_rk4 = np.arange(0, 2 * T, dt)        # tiempo de evaluacion rk4
y_rk4 = np.zeros((len(t_eval_rk4), 2))

# Cond Inicial
y_rk4[0] = [phi_0, phipunto_0]

# función que rellena la matriz y_rk4 con los valores obtenidos mediante rk4
for i in range(1, len(t_eval_rk4)):
    y_rk4[i] = paso_rk4(pendulo_viscoso, dt, t_eval_rk4[i-1], y_rk4[i-1])


plt.figure(1)
plt.clf()

plt.plot(t_to_plot, phi_pequenas_osc, label='Pequeñas osc.')
plt.plot(t_eval_rk4, y_rk4[:, 0], label='Solución RK4')

plt.xlabel('Tiempo [s]', fontsize=15)
plt.ylabel(r'$\phi(t)$', fontsize=15)
plt.legend()
plt.show()
