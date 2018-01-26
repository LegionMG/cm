import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


import scipy.sparse
import scipy.sparse.linalg


def func(x, y, t):
    return 2 * np.sin(x + y) + (2 * t) / (t ** 2 + 1)


def B(x, t):
    return np.sin(x) + np.log(t ** 2 + 1)
    
    
def I(x, y):
    return np.sin(x + y)




def R(x, y, t):
    return np.sin(x + y) + np.log(t ** 2 + 1)


def exact(Nx, Nt,L, T):
    
    x   = np.linspace(0, L, Nx+1)
    y   = np.linspace(0, L, Nx+1) 
    t   = np.linspace(0, T, Nt+1)

    U_real = np.zeros(shape=(Nx + 1, Nx + 1, Nt + 1))
    for k in range(Nt + 1):
        for j in range(Nx + 1):
            for i in range(Nx + 1):
                U_real[i, j, k] = R(x[i], y[j], t[k])
    return x, y, t, U_real


def error_max(U, h, delta, Nx, Nt, L, T):
    """
    Вычисление ошибки
    """
    x, y, t, U_real = exact(Nx, Nt, L, T)
    return np.max(np.abs(U - U_real))

def plot_t(x, y, t, U):    
    """
    Построение графика U(x, y) при t 
    """
    title = "U(x, y, t=" + str(t) + ")"
    fig = plt.figure()
    ax = fig.gca(projection="3d")
    X, Y = np.meshgrid(x, y)
    Z = U[:, :, t]
    ax.plot_surface(X, Y, Z, cmap="viridis", rstride=1, cstride=1, linewidth=0)    
    ax.view_init(30, 220);
    plt.suptitle(title)
    plt.show()
    
def plot_x(y, t, x, U):
    """
    Построение графика U(y, t) при x 
    """
    title = "U(y, t, x=" + str(x) + ")"
    fig = plt.figure()
    ax = fig.gca(projection="3d")
    X, Y = np.meshgrid(y, t)
    Z = U[x, :, :].T
    ax.plot_surface(X, Y, Z, cmap="viridis", rstride=1, cstride=1, linewidth=0)    
    ax.view_init(30, 220);
    plt.suptitle(title)
    plt.show() 
       
def result(schema):
    _errors = []
    _steps = []
    _deltas = []
    steps = [np.pi / 5, np.pi / 10, np.pi / 20]
    deltas = [0.25, 0.10, 0.05]

    a   = 1
    L   = 2*np.pi
    T   = 10

    Nxs = [5, 10, 20]
    Nts = [20, 40, 80]

    dxs = [L/x for x in Nxs]
    dts = [T/x for x in Nts]
    
    x, y, t, U, error = schema(dxs[0], dts[0], Nxs[0], Nts[0], L, T)
    _, _, _, U_real = exact(Nxs[0], Nts[0], L, T)
    
    print("Exact:")
    plot_t(x, y, -1, U_real)
    print("Alt_dim:")
    plot_t(x, y, -1, U)
    
    for i, h in enumerate(dxs):
        for j, delta in enumerate(dts):
            x, y, t, U, err = schema(h, delta, Nxs[i], Nts[j], L, T)
            _errors.append(err)
            _steps.append(h)
            _deltas.append(delta)

    df = pd.DataFrame(data={
        "step": _steps, 
        "delta": _deltas,
        "error": _errors
    })
    
    print("Max errors:")
    return df

def alt_dim(dx, dt, Nx, Nt, L, T):
    dt = dt / 2
    Nt *= 2

    F = dt / (dx ** 2)

    x   = np.linspace(0, L, Nx+1)
    y   = np.linspace(0, L, Nx+1) 
    t   = np.linspace(0, T, Nt+1)

    u_s = np.zeros(shape=(Nx + 1, Nx + 1, Nt + 1))

    for j in range(Nx + 1):
        u_s[j,:,0] = I(x[j], y)

    #не вышло векторизовать :()
    for i in range(Nx + 1):
        u_s[0, i, :] = B(y[i], t)
        u_s[-1, i, :] = B(y[i], t)
        u_s[i, 0, :] = B(x[i], t)
        u_s[i, -1, :] = B(x[i], t)



    b = np.zeros(shape=(Nx + 1, Nx + 1))

    main  = np.zeros(Nx+1)
    lower = np.zeros(Nx)
    upper = np.zeros(Nx)

    main[:]  = 1 + 2 * F
    lower[:] = -F
    upper[:] = -F
    upper[0] = 0
    lower[-1] = 0
    main[0] = main[-1] = 1

    A = scipy.sparse.diags(
            diagonals=[main, lower, upper],
            offsets=[0, -1, 1], shape=(Nx+1, Nx+1),
            format='csr')

    for i in range(Nt):
        if i % 2 != 0:
            # Прогонка по x

            for j in range(1, Nx):
                b[:, j] = u_s[:,j, i] + F * (u_s[:, j-1, i] - 2*u_s[:, j, i] + u_s[:, j+1, i]) + dt*func(x, y[j], t[i+1])
            for j in range(Nx):
                u_s[j, :, i + 1] = scipy.sparse.linalg.spsolve(A, b[j])
         
        else:
            # Прогонка по y
            for j in range(1, Nx):
                b[j] = u_s[j, :, i] + F*(u_s[j-1, :, i] - 2*u_s[j, :, i] + u_s[j+1, :, i]) + dt*func(x[j], y, t[i+1])
            for j in range(Nx):
                u_s[:, j, i + 1] = scipy.sparse.linalg.spsolve(A, b[j])   
        
        #Должно быть не нужно, всё должно вычисляться через b
        u_s[:, 0, i + 1] = B(x, t[i + 1])
        u_s[:, -1, i + 1] = B(x, t[i + 1])
        u_s[0, :, i + 1] = B(y, t[i + 1])
        u_s[-1, :, i + 1] = B(y, t[i + 1])

    real = np.zeros(shape=(Nx + 1, Nx + 1, Nt + 1))
    for j in range(Nx + 1):
        for i in range(Nx + 1):
            real[i, j, :] = R(x[i], y[j], t)


    #print(R(x,y,t[1]))
    error = np.max(np.abs(u_s - real))

    return x, y, t[::2], u_s[:, :, ::2], error







df = result(alt_dim)

print(df)