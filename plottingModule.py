import fieldsModule
from constants import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # This import is required for the kwarg projection='3d' in the call to fig.add_subplot


def plot3dWaist(index):
    z = 1.e-5
    t = z/c
    xs = np.array([xmin+xi*dx for xi in range(grid+1)])/w0
    ys = np.array([ymin+yi*dy for yi in range(grid+1)])/w0
    xsMat, ysMat = np.meshgrid(xs, ys)  # turns xs and ys into matrices: loop over xs inside of loop over ys
    zsMat = np.array([[0.0 for ix in range(grid+1)] for iy in range(grid+1)])
    for iy in range(grid+1):
        y = ymin+iy*dy
        for ix in range(grid+1):
            x = xmin+ix*dx
            fieldsModule.makeFields(x, y, z, t)
            zsMat[iy][ix] = abs(fieldsModule.ceb[index])
    #
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(xsMat, ysMat, zsMat)
    ax.set_xlabel('x/w0')
    ax.set_ylabel('y/w0')
    ax.set_zlabel('abs(ceb['+str(index)+'])')
    plt.show()
    return


def plot2dWaist(index):
    z = 1.e-5
    t = z/c
    xs = np.array([xmin+xi*dx for xi in range(grid+1)])/w0
    ys = np.array([ymin+yi*dy for yi in range(grid+1)])/w0
    xsMat, ysMat = np.meshgrid(xs, ys)  # turns xs and ys into matrices: loop over xs inside of loop over ys
    zsMat = np.array([[0.0 for ix in range(grid+1)] for iy in range(grid+1)])
    for iy in range(grid+1):
        y = ymin+iy*dy
        for ix in range(grid+1):
            x = xmin+ix*dx
            fieldsModule.makeFields(x, y, z, t)
            zsMat[iy][ix] = abs(fieldsModule.ceb[index])
    #
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.contourf(xsMat, ysMat, zsMat, 25)
    plt.show()
    return


def plot1dWaist(index):
    z = 1.e-5
    t = z/c
    xs = np.array([xmin+xi*dx for xi in range(grid+1)])/w0
    zs = np.array([0.0 for ix in range(grid+1)])
    for ix in range(grid+1):
        x = xmin+ix*dx
        fieldsModule.makeFields(x, 1.e-5, z, t)
        zs[ix] = abs(fieldsModule.ceb[index])
    #
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(xs, zs)
    plt.show()
    return
