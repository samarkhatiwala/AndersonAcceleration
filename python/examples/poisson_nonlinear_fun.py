import numpy as np
from scipy.sparse import spdiags, kron, eye


def poisson_nonlinear_fun(f, bx0, bxf, by0, byf, D, nx, ny):

    x0, xf, y0, yf = tuple(D)

    hx = (xf - x0) / (nx + 1)
    hy = (yf - y0) / (ny + 1)

    x = np.linspace(x0, xf, nx + 2)
    y = np.linspace(y0, yf, ny + 2)

    X, Y = np.meshgrid(x, y)
    X = np.transpose(X)
    Y = np.transpose(Y)

    Iint = np.arange(1, nx + 1)
    Jint = np.arange(1, ny + 1)
    Xint = X[Iint, :][:, Jint]
    Yint = Y[Iint, :][:, Jint]

    rhs = f(Xint, Yint)
    rhs[:, 0] += by0(X[Iint, 0]) / (hy ** 2)
    rhs[:, ny - 1] += byf(X[Iint, -1]) / (hy ** 2)
    rhs[0, :] += bx0(Y[0, Jint]) / (hx ** 2)
    rhs[nx - 1, :] += bxf(Y[-1, Jint]) / (hx ** 2)

    b = np.reshape(rhs, (nx * ny, 1), order="F")

    Ix = eye(nx)
    Iy = eye(ny)
    ex = np.ones((nx, 1))
    ey = np.ones((ny, 1))

    data =np.array([ex/(hx**2), -2*ex*(1/hx**2+1/1/hy**2), ex/(hx**2)])
    data = np.squeeze(data)
    T = spdiags(data, [-1, 0, 1], m=nx, n=nx)
    S = spdiags(np.squeeze(np.array([ey, ey])), [-1, 1], m=ny, n=ny)

    A = -(kron(Iy, T) + kron(S, Ix/(hy**2)))

    return A, b