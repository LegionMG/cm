from numpy import *
import matplotlib.pyplot as plt
from scipy.linalg import norm

def initFun(tInit, x):
    """This function calculates initial function and is used to fulfill
    the first (presely, zero) time layer in the grid solution aray
    """
    res = sin(x) * cos(tInit)
    return res

def bFun(t, x):
    """This function calculates left boundary function
    """
    res = sin(xLeft) * cos(t)
    return res


def nFun(u):
    """This function calculates nonlinear part of difference operator
    """
    return exp(u)

def nFunDer(u):
    """This function calculates the derivative with respect to time for
    nonlinear part of difference operator
    """
    return exp(u)

def heterFun(t, x):
    """This function calculates the heterogenious function
    """
    res = -exp(sin(x)*cos(t))*sin(x)*sin(t)+sin(x)*cos(t)
    return res

def realS(t, x):
    """This function calculates the exact solution and
    is used to estimate the error
    """
    res = sin(x)*cos(t)
    return res

xLeft = 0            #xLeft - left boundary 
xRight = 3.14159265  #xRight - right boundary, we solve the boundary problen on the segment [xLeft, xRight]
tInit = 0            #tInit - initial time
tFinal = 8.0         #tFinal - final time, we solve the boundary problen on the segment [tInint, tFinal]
Nx = 31  #numSeg - number of subsegments in segment [x_left, x_right]
Nt = 80  #numSeg - number of subsegments in segment [tInint, tFinal]
errorNorm = 0
epsNewton = 10**(-5)

#Below we calculate the time step dt,
#construct the uniform time grid and allocate memory for solution-array
#Note, that number of dots in time grid = numSeg + 1
dx = (xRight - xLeft)/Nx
sG = linspace(xLeft, xRight, Nx + 1)

dt = (tFinal - tInit)/Nt
tG = linspace(tInit, tFinal, Nt + 1)

uS = zeros((Nt + 1, Nx + 1))
uError = zeros((Nt + 1, Nx + 1))
alpha = zeros(Nx + 1)
beta = zeros(Nx + 1)
yk = zeros(Nx + 1)
yk1= zeros(Nx + 1)

F = dt/(2*dx**2)

for i in range(Nx + 1):
    uS[0][i] = initFun(tInit, sG[i])


for j in range(Nt):
    uS[j+1][0] = bFun(tG[j+1], xLeft)
    uS[j+1][Nx] = bFun(tG[j+1], xRight)
    flag = True
    yk = copy(uS[j])
    alpha[0] = 0
    beta[0] = bFun(tG[j+1], xLeft)
    yk1[Nx] = bFun(tG[j+1], xRight)
    numOfIterations = 0
    while flag:
        numOfIterations += 1
        for i in range(0, Nx):
            Fi = -(nFun(yk[i]) - nFun(uS[j][i]) - nFunDer(yk[i])*yk[i] - dt*heterFun(tG[j], sG[i]) - dt/(2*dx**2) * (uS[j][i-1] - 2*uS[j][i] + uS[j][i+1]))
            Ci = nFunDer(yk[i]) + dt/dx**2
            alpha[i+1] = F/(Ci - alpha[i]*F)
            beta[i+1] = (F*beta[i] + Fi)/(Ci - alpha[i]*F)
        for i in range(Nx -1, 0, -1 ):
            yk1[i] = alpha[i+1]*yk1[i+1] + beta[i+1]
        flag = (norm(yk1 - yk) > epsNewton)
        yk = copy(yk1)
    uS[j+1] = copy(yk1)

    for i in range(0, Nx + 1):
        uError[j+1][i] = abs(uS[j+1][i] - realS(tG[j+1], sG[i]))
        if (uError[j+1][i] > errorNorm):
            errorNorm = uError[j+1][i]

#print(uError)
#print(uS)
print(errorNorm)
print(numOfIterations)
#plt.plot(tG, uError[:,16])
#print(uS[:, -1])
plt.plot(tG, uS[:, -2])

plt.show()