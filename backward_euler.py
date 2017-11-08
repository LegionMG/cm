import numpy as np
import matplotlib.pyplot as plt
import math
from matplotlib.widgets import Button

import scipy.sparse
import scipy.sparse.linalg

class Index(object):
    ind = 0

    def next(self, event):
        self.ind = (self.ind + 1) % len(u_s)
        s = u_s[self.ind]
        s_b = u_s[self.ind]
        ax.cla()
        ax.plot(t, s_b[0], 'bo')
        ax.plot(t, s[1], 'r--')
        plt.draw()


    def prev(self, event):
        self.ind = (self.ind - 1) % len(u_s)
        s = u_s[self.ind]
        s_b = u_s[self.ind]
        ax.cla()
        ax.plot(t, s_b[0], 'bo')
        ax.plot(t, s[1], 'r--')
        plt.draw()


def I(x):
    return math.sin(x)

def G(x, t):
    return np.sin(x) + 2*t/(t**2 + 1)

def B(t):
    return math.log(t**2 + 1)

a   = 1
L   = math.pi
T   = 10
Nx  = 5
Nt  = 40

x   = np.linspace(0, L, Nx+1) 
dx  = x[1] - x[0]
t   = np.linspace(0, T, Nt+1)
dt  = t[1] - t[0]
F   = a*dt/dx**2
u   = np.zeros(Nx+1)
u_1 = np.zeros(Nx+1)


main  = np.zeros(Nx+1)
lower = np.zeros(Nx)
upper = np.zeros(Nx)
b     = np.sin(x)


main[:]  = 1 + 2*F
lower[:] = -F 
upper[:] = -F 


u_s = []


#решение
for n in t:
    main[0] = main[Nx] = B(n)
    A = scipy.sparse.diags(
        diagonals=[main, lower, upper],
        offsets=[0, -1, 1], shape=(Nx+1, Nx+1),
        format='csr')
    u_s.append((np.copy(u_1), np.sin(x) + math.log(n**2 + 1)))
    b = u_1
    b[0] = b[-1] = B(n) 
    print(b)
    u[:] = scipy.sparse.linalg.spsolve(A, b)
    u_1[:] = u + dt * (np.sin(x) + 2*n/(n**2 + 1)) 


#Отрисовка
fig, ax = plt.subplots()
callback = Index()

axprev = plt.axes([0.7, 0.05, 0.1, 0.075])
axnext = plt.axes([0.81, 0.05, 0.1, 0.075])
bnext = Button(axprev, 'Next')
bnext.on_clicked(callback.next)
bprev = Button(axnext, 'Previous')
bprev.on_clicked(callback.prev)

plt.subplots_adjust(bottom=0.2)
t = x
s = u_s[0]
s_b = u_s[-1]
our, = ax.plot(t, s_b[0], 'bo')
real, = ax.plot(t, s[1], 'r--')



plt.show()