"""

"""

import numpy as np 
import matplotlib.pyplot as plt 

g = 9.8 # m/s^2
l = 5.533 # m
gamma = 2.533 # s^-1


phi_0 = 1 
omega_0 = 1

def pendulo_viscoso(t,y):
    output = np.array(y[1], -gamma * y[1] - (g/l)* np.sin(y[0]))
    return output








































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