import numpy as np
import matplotlib.pyplot as plt
 
# constantes
p_phi = 10
g = 9.8
m = 1
l = 1
 
theta_to_plot = np.linspace(-4*np.pi , 4*np.pi , 10000000)

def func(x):
  u_eff = ((p_phi**2) / (2*m*l*l*np.sin(theta_to_plot)*np.sin(theta_to_plot))) + m*g*l*np.cos(theta_to_plot)
  return u_eff


plt.figure(1)
plt.clf()
plt.plot(theta_to_plot, func(theta_to_plot))
 
plt.xlabel(r'$ \theta $', fontsize=15)
plt.ylabel(r'$U_{eff}(\theta)$', fontsize=15)
plt.show()