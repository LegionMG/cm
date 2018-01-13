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
print(A.toarray())

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

print(errors)
t = x
s = u_s[0]
s_b = u_s[0]
our, = ax.plot(t, s_b[0], 'bo')
real, = ax.plot(t, s[1], 'r--')

plt.show()