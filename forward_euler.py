
import matplotlib.pyplot as plt
from matplotlib.widgets import Button
from config import *

class Index(object):
    ind = 0

    def next(self, event):
        self.ind = (self.ind + 1) % len(u_s)

        s = u_s[self.ind]
        ax.cla()
        ax.plot(t, s[0], 'bo')
        ax.plot(t, s[1], 'r--')
        plt.draw()
    def prev(self, event):
        self.ind = (self.ind - 1) % len(u_s)
        
        s = u_s[self.ind]
        ax.cla()
        ax.plot(t, s[0], 'bo')
        ax.plot(t, s[1], 'r--')
        plt.draw()



Nx  = int(L/dx)
Nt  = int(T/dt)
x   = np.linspace(0, L, Nx+1) 
t   = np.linspace(0, T, Nt+1)
F   = a*dt/dx**2
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
our, = ax.plot(t, s[0], 'bo')
real, = ax.plot(t, s[1], 'r--')

print(errors)
plt.show()