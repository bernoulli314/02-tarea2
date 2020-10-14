"""

"""

import numpy as np 
import matplotlib.pyplot as plt 

g = 9.8 # m/s^2
l = 5.533 # m
gamma = 2.533 # s^-1


phi_0 = np.pi / gamma
phipunto_0 = 0

# Graficamos la solución de pequeñas oscilaciones
# La sol. de p.o. es phi(t) = A * e^(-1/2 * t * gamma) * cos(omega * t + phi_0)
# Con A} cte dependiente de phi_0

omega = np.sqrt(g/l - (gamma **2) / 4)
A = phi_0 / np.cos(phi_0)
T = 2 * np.pi / omega

t_to_plot = np.linspace(0 , 2 * T, 200)

phi_pequenas_osc = A * np.e ** (- t_to_plot * gamma / 2) * np.cos(omega * t_to_plot + phi_0)

plt.figure(1)
plt.clf()

plt.plot(t_to_plot, phi_pequenas_osc, label='Pequ. osc.')

plt.xlabel('Tiempo[s]')
plt.ylabel(r'$\phi(t)$', fontsize=15)
plt.legend()
plt.show()



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

T = 2*np.pi / omega

dt = 0.001
t_eval_rk4 = np.arange(0, 2 * T, dt)
y_rk4 = np.zeros((len(t_eval_rk4), 2))

# cond inicial
y_rk4[0] = [phi_0, phipunto_0]
for i in range(1, len(t_eval_rk4)):
    y_rk4[i] = paso_rk4(pendulo_viscoso, dt, t_eval_rk4[i-1], y_rk4[i-1])



plt.plot(t_eval_rk4, y_rk4[:,0], label='Sol. RK4')

plt.xlabel('Tiempo [s]', fontsize=15)
plt.ylabel(r'$\phi(t)$', fontsize=15)
plt.legend()
plt.show()