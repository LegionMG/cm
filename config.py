import numpy as np
import math

def I(x):
    return math.sin(x)

def G(x, t):
    return np.sin(x) + 2*t/(t**2 + 1)

def B(t):
    return math.log(t**2 + 1)

a   = 1
L   = math.pi
T   = 10

dx  = L/5
dt  = 1/40
Nx  = int(L/dx)
Nt  = int(T/dt)


dxs = [L/5, L/10, L/20]
dts = [1/10,1/40, 1/70]



#Настройка
make_graphs = True #два режима - рисовать график\прогонять на нескольких вариантах шагов. True\False соответственно.
graph = "bw" #fw, bw, kn(forward euler, backward euler, krank-nicholson)