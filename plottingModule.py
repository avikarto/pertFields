import numpy as np
import fieldsModule
import constants as const
import matplotlib.pyplot as plt
# This import is required for the kwarg projection='3d' in the call to fig.add_subplot
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

xmin, dx = const.xmin, const.dx
ymin, dy = const.ymin, const.dy


def plot3dWaist(field_type_index):
	z = 1.e-5
	t = z / const.c
	xs = np.array([xmin + (xi * dx) for xi in range(const.grid + 1)]) / const.w0
	ys = np.array([ymin + (yi * dy) for yi in range(const.grid + 1)]) / const.w0
	# Turn xs and ys into matrices; loop over xs inside of loop over ys
	xsMat, ysMat = np.meshgrid(xs, ys)
	zsMat = np.array([[0.0 for ix in range(const.grid + 1)] for iy in range(const.grid + 1)])
	for iy in range(const.grid + 1):
		y = ymin + (iy * dy)
		for ix in range(const.grid + 1):
			x = xmin + (ix * dx)
			ceb = fieldsModule.makeFields(x, y, z, t)
			zsMat[iy][ix] = abs(ceb[field_type_index])
	# for iy

	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	ax.plot_surface(xsMat, ysMat, zsMat)
	ax.set_xlabel('x/w0')
	ax.set_ylabel('y/w0')
	ax.set_zlabel('abs(ceb[' + str(field_type_index) + '])')
	plt.show()
	return


def plot2dWaist(field_type_index):
	z = 1.e-5
	t = z / const.c
	xs = np.array([xmin + (xi * dx) for xi in range(const.grid + 1)]) / const.w0
	ys = np.array([ymin + (yi * dy) for yi in range(const.grid + 1)]) / const.w0
	# Turn xs and ys into matrices; loop over xs inside of loop over ys
	xsMat, ysMat = np.meshgrid(xs, ys)
	zsMat = np.array([[0.0 for ix in range(const.grid + 1)] for iy in range(const.grid + 1)])
	for iy in range(const.grid + 1):
		y = ymin + (iy * dy)
		for ix in range(const.grid + 1):
			x = xmin + (ix * dx)
			ceb = fieldsModule.makeFields(x, y, z, t)
			zsMat[iy][ix] = abs(ceb[field_type_index])
	# for iy

	fig = plt.figure()
	ax = fig.add_subplot(111)
	ax.contourf(xsMat, ysMat, zsMat, 25)
	plt.show()
	return


def plot1dWaist(field_type_index):
	z = 1.e-5
	t = z / const.c
	xs = np.array([xmin + (xi * dx) for xi in range(const.grid + 1)]) / const.w0
	zs = np.array([0.0 for ix in range(const.grid + 1)])
	for ix in range(const.grid + 1):
		x = xmin + (ix * dx)
		ceb = fieldsModule.makeFields(x, 1.e-5, z, t)
		zs[ix] = abs(ceb[field_type_index])
	# for ix

	fig = plt.figure()
	ax = fig.add_subplot(111)
	ax.plot(xs, zs)
	plt.show()
	return
