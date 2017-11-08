import numpy as np
import matplotlib.pyplot as plt
import math

def I(x):
    return math.sin(x)
def G(x, t):
    return math.sin(x) + 2*t/(t**2 + 1)

a = 1
L = math.pi
T = 10
Nx = 5
Nt = 40

x = np.linspace(0, L, Nx+1)    # mesh points in space
dx = math.pi/5
t = np.linspace(0, T, Nt+1)    # mesh points in time
print(x)
dt = 0.25
F = a*dt/dx**2
u   = np.zeros(Nx+1)           # unknown u at new time level
u_1 = np.zeros(Nx+1)           # u at the previous time level

for i in range(0, Nx+1):
    u_1[i] = I(x[i])

for n in t:
    plt.figure(1)
    plt.plot(x, u_1, 'bo')
    plt.plot(x, np.sin(x) + math.log(n**2 + 1), 'r--') # real function
    plt.show()

    # Compute u at inner mesh points
    for i in range(1, Nx):
        u[i] = u_1[i] + F*(u_1[i-1] - 2*u_1[i] + u_1[i+1])+dt*G(x[i], n)

    # Insert boundary conditions
    u[0] = math.log((n)**2 + 1);  u[Nx] = math.log((n)**2 + 1)

    # Update u_1 before next step
    u_1[:]= u

