import matplotlib.pyplot as plt
from matplotlib.widgets import Button

from config import *

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


Nx  = int(L/dx)
Nt  = int(T/dt)
x   = np.linspace(0, L, Nx+1) 
t   = np.linspace(0, T, Nt+1)
F   = a*dt/dx**2


def krank_nicholson():
    u   = np.zeros(Nx+1)
    u_1 = np.sin(x)

    main  = np.zeros(Nx+1)
    lower = np.zeros(Nx)
    upper = np.zeros(Nx)
    b     = np.sin(x)


    main[:]  = 1 + 2*F
    lower[:] = -F 
    upper[:] = -F 
    A = scipy.sparse.diags(
        diagonals=[main, lower, upper],
        offsets=[0, -1, 1], shape=(Nx+1, Nx+1),
        format='csr')

    u_s = []
    u_b = np.zeros(Nx+1)
    u_f = np.zeros(Nx+1)

    errors = []
    #решение
    for n in t:
        z = np.sin(x) + math.log(n**2 + 1)
        errors.append(max(np.abs(u_1 - z)))
        if n != t[0]:
            b = (u_1 + dt * (np.sin(x) + 2*n/(n**2 + 1))) 
        b[0] = b[Nx] = B(n) 
        u[:] = scipy.sparse.linalg.spsolve(A, b)
        u_b[:] = u


        for i in range(1, Nx):
            u_f[i] = u_1[i] + F*(u_1[i-1] - 2*u_1[i] + u_1[i+1])+dt*G(x[i], n)

        u_f[0] = math.log((n)**2 + 1);  u_f[Nx] = math.log((n)**2 + 1)

        u_1[:] = 0.5*u_f + 0.5*(u_b + dt * (np.sin(x) + 2*n/(n**2 + 1))) 


        u_s.append((np.copy(u_1), np.sin(x) + math.log(n**2 + 1)))
       # print(u_f, u_b, u_1, np.sin(x) + math.log(n**2 + 1))

    return u_s, max(errors)

def forward_euler():
    u   = np.zeros(Nx+1)           
    u_1 = np.sin(x)         
    u_s = []


    errors = []

    for n in t:
        z = np.sin(x) + math.log(n**2 + 1)
        errors.append(max(np.abs(u_1 - z)))
        u_s.append((np.copy(u_1), np.sin(x) + math.log(n**2 + 1)))
        for i in range(1, Nx):
            u[i] = u_1[i] + F*(u_1[i-1] - 2*u_1[i] + u_1[i+1])+dt*G(x[i], n)

        u[0] = math.log((n)**2 + 1);  u[Nx] = math.log((n)**2 + 1)

        u_1[:] = u

    return u_s, max(errors)


def backward_euler():
    u   = np.zeros(Nx+1)
    u_1 = np.sin(x)


    main  = np.zeros(Nx+1)
    lower = np.zeros(Nx)
    upper = np.zeros(Nx)
    b     = np.sin(x)


    main[:]  = 1 + 2*F
    lower[:] = -F 
    upper[:] = -F 
    A = scipy.sparse.diags(
        diagonals=[main, lower, upper],
        offsets=[0, -1, 1], shape=(Nx+1, Nx+1),
        format='csr')

    u_s = []
    errors = []

    #решение
    for n in t:
        z = np.sin(x) + math.log(n**2 + 1)
        errors.append(max(np.abs(u_1 - z)))
        if n != t[0]:
            b = (u_1 + dt * (np.sin(x) + 2*n/(n**2 + 1))) 
        b[0] = b[Nx] = B(n) 
        u[:] = scipy.sparse.linalg.spsolve(A, b)
        u_1[:] = u
        
        u_s.append((np.copy(u_1), np.sin(x) + math.log(n**2 + 1)))


    return u_s, max(errors)









#Отрисовка
if make_graphs:
    if graph == 'fw':
        u_s, e = forward_euler()
    elif graph == 'bw':
        u_s, e = backward_euler()
    else:
        u_s, e = krank_nicholson()
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
    s = u_s[0]
    s_b = u_s[0]
    our, = ax.plot(t, s_b[0], 'bo')
    real, = ax.plot(t, s[1], 'r--')
    plt.show()
else:
    data = []
    for dx in dxs:
        for dt in dts:
            Nx  = int(L/dx)
            Nt  = int(T/dt)
            x   = np.linspace(0, L, Nx+1) 
            t   = np.linspace(0, T, Nt+1)
            F   = a*dt/dx**2

            u_s, e = forward_euler()

            data.append((dt, dx, e))

    print("forward_euler:")
    for line in data:
        print(' '.join(map('{0:.4f}'.format, line)))
    print("#"*20)

    data = []
    for dx in dxs:
        for dt in dts:
            Nx  = int(L/dx)
            Nt  = int(T/dt)
            x   = np.linspace(0, L, Nx+1) 
            t   = np.linspace(0, T, Nt+1)
            F   = a*dt/dx**2

            u_s, e = backward_euler()

            data.append((dt, dx, e))

    print("backward_euler:")
    for line in data:
        print(' '.join(map('{0:.4f}'.format, line)))
    print("#"*20)



    data = []
    for dx in dxs:
        for dt in dts:
            Nx  = int(L/dx)
            Nt  = int(T/dt)
            x   = np.linspace(0, L, Nx+1) 
            t   = np.linspace(0, T, Nt+1)
            F   = a*dt/dx**2

            u_s, e = krank_nicholson()

            data.append((dt, dx, e))

    print("krank_nicholson:")
    for line in data:
        print(' '.join(map('{0:.4f}'.format, line)))
    print("#"*20)

