import numpy as np
import math

def I(x):
    return np.sin(x)

def G(x, t):
    return np.sin(x) + (2*t)/(t**2 + 1)

def B(t):
    return np.log(t**2 + 1)

def R(x, t):
    return I(x) + B(t)

a   = 1
L   = math.pi
T   = 10



Nx  = 5
Nt  = 40

dx  = L/Nx
dt  = T/Nt


Nxs = [5, 10, 20, 40, 70]
Nts = [40, 80, 150, 300, 600]

dxs = [L/x for x in Nxs]
dts = [T/x for x in Nts]


#Настройка
make_graphs = False #два режима - рисовать график\прогонять на нескольких вариантах шагов. True\False соответственно.
graph = "kn" #fw, bw, kn(forward euler, backward euler, krank-nicholson)