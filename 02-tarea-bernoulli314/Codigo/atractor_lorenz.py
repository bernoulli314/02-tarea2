"""

"""

import numpy as np 
import matplotlib.pyplot as plt 
from scipy.integrate import RK45
from mpl_toolkits.mplot3d import Axes3D

sigma = 10
beta = 8/3
rho = 28

def lorenz(s,w):
    output = np.array([sigma * (w[1] - w[0]) , w[0] * (rho - w[2]) - w[1] , w[0] * w[1] - beta * w[2]])
    return output

w_0 = np.array([1,1,1])

rk4 = RK45(lorenz, 0.0, w_0, 0.001)