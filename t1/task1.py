import matplotlib.pyplot as plt
from matplotlib.widgets import Button

from config import *

import scipy.sparse
import scipy.sparse.linalg

class Index(object):
    ind = 0

    def next(self, event):
        self.ind = (self.ind + 1) % len(meths)
        s = reals[self.ind]
        s_b = meths[self.ind]
        ax.cla()
        ax.plot(t, s_b, 'bo')
        ax.plot(t, s, 'r--')
        plt.draw()


    def prev(self, event):
        self.ind = (self.ind - 1) % len(meths)
        s = reals[self.ind]
        s_b = meths[self.ind]
        ax.cla()
        ax.plot(t, s_b, 'bo')
        ax.plot(t, s, 'r--')
        plt.draw()


Nx  = int(L/dx)
Nt  = int(T/dt)
x   = np.linspace(0, L, Nx+1) 
t   = np.linspace(0, T, Nt+1)
F   = a*dt/dx**2


def theta_met(theta):
    u   = np.zeros(Nx+1)
    u_1 = np.sin(x)

    main  = np.zeros(Nx+1)
    lower = np.zeros(Nx)
    upper = np.zeros(Nx)
    b     = np.sin(x)

    #theta = 0.5

    main[:]  = 1 + 2*F *theta
    lower[:] = -F * theta
    upper[:] = -F * theta
    A = scipy.sparse.diags(
        diagonals=[main, lower, upper],
        offsets=[0, -1, 1], shape=(Nx+1, Nx+1),
        format='csr')

    #print(u_s)
    meths = [[] for _ in range(len(t))]
    reals = [(np.sin(x) + math.log(i**2 + 1))[:] for i in t]
    errors = []
    #решение
    for i in range(len(t)):
        

        b[1:-1] = u_1[1:-1] + F * (1-theta) *(u_1[:-2] - 2*u_1[1:-1] + u_1[2:]) + dt*G(x[1:-1], t[i])
        b[0] = b[-1] = B(t[i]) 
        u[:] = scipy.sparse.linalg.spsolve(A, b)
        u_1[:] = u


        meths[i] = (np.copy(u_1))
        z = np.sin(x) + math.log(t[i]**2 + 1)
        errors.append(max(np.abs(u_1 - z)))

    #print(u_s)
    return reals, meths, max(errors)


d = {"fw": 0, "bw": 1, "kn": 0.5}


#Отрисовка
if make_graphs:
    reals, meths, e = theta_met(d[graph])
    #print(meths)
    print("Error: {0}".format(e))
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
    print(meths[0], reals[0])
    our, = ax.plot(t, list(meths[0]), 'bo')
    real, = ax.plot(t, reals[0], 'r--')
    plt.show()
else:
    for k in d.keys():
        data = []
        for i in range(len(dxs)):
            dt = dts[i]
            dx = dxs[i]
            Nx  = int(L/dx)
            Nt  = int(T/dt)
            x   = np.linspace(0, L, Nx+1) 
            t   = np.linspace(0, T, Nt+1)
            F   = a*dt/dx**2

            reals, meth, e = theta_met(d[k])

            data.append((dt, dx, e))
        print(k)
        for line in data:
            print(' '.join(map('{0:.4f}'.format, line)))
        print("#"*20)

    
