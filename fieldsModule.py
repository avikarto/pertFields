from math import pi, sqrt
from constants import (
	nn, mm, nnReal, mmReal,
	ii, c, c2,
	phi0, omega0, s, sReal, pertOrder, zr, epsilonc2, polar
)
from mathModule import binomial, fact, my_gamma
from numpy import arctan2, exp, array, nan


def lambdaNM(phi=0.0, xi=0.0j, beta=0.0j):  # \Lambda_{n,m}
	res = (
		(-1.0)**(nn + mm)
		* 2.0**((2 * nn) + mm)
		* sqrt(2.0 * pi)
		* fact(nn)
		* exp(ii * phi0)
		* xi**(mm / 2.0)
		* beta**(-nn - (mm / 2.0) - 1)
		* exp(ii * mm * phi)
	)
	return res


def kappa(alpha, delta):  # \kappa_{\alpha,\delta}
	return (-1.0)**(delta - alpha) * binomial(2 * alpha, 2 * alpha - delta) / fact(delta - alpha)


def bigG(v1, v2, v3):  # G_{n,m,j}
	return (-1.0)**v3 * fact(v1 + v2) / (fact(v1 - v3) * fact(v2 + v3) * fact(v3))


def cAD(alpha, delta, jay):  # c_{\alpha,\delta}
	# Generate a single cAD value
	gamma = mmReal / 2.0 + sReal + float(jay)
	res = (
		kappa(alpha, delta)
		* bigG(nn + delta, mm, jay)
		* (fact(nn + delta) / fact(nn))
		* omega0**alpha
		* (sReal / omega0)**(sReal - gamma + float(alpha))
		* my_gamma(gamma + 1 - alpha) / my_gamma(s + 1)
	)
	return res


def make_cADs_matrix():
	# Generating a matrix of cAD values
	# Reference as cADs[alpha][delta][j]
	cADs = [
		[
			[nan for j in range(0, nn + 2 * pertOrder + 1)]
			for d in range(0, 2 * pertOrder + 1)
		]
		for a in range(0, pertOrder + 1)
	]
	for a in range(0, pertOrder + 1):
		for d in range(a, 2 * a + 1):
			for j in range(0, nn + d + 1):
				cADs[a][d][j] = cAD(a, d, j)
	return cADs


def makeFields(x, y, z, t) -> list:
	''' Calculate the phasor and complex E/B fields at any space-time coordinate.
	-----
	Accepts: x, y, z, t - the space-time coordiantes.
	Returns: ceb - complex values for [phasor,E1,E2,E3,B1,B2,B3]
	'''

	rho = sqrt(x**2 + y**2)
	if rho < 1.e-5:
		rho = 1.e-5
	phi = arctan2(y, x)  # radians
	beta = 1.0 + ii * z / zr
	xi = rho**2 / (2.0 * c * beta * zr)
	bigT = 1.0 + omega0 * ((-ii * z / c) + xi + (ii * t)) / sReal
	lambda_nm = lambdaNM(phi=phi, xi=xi, beta=beta)

	# calculating the time-domain phasor, and its derivatives
	ceb = array([0.0j for i in range(7)])
	dUdT = 0.0j
	d2UdT2 = 0.0j
	dUdxi = 0.0j
	d2Udxi2 = 0.0j
	dUdb = 0.0j
	d2UdxidT = 0.0j
	d2Udxidb = 0.0j
	d2UdbdT = 0.0j

	cADs = make_cADs_matrix()

	for iAlpha in range(0, pertOrder + 1):
		for iDelta in range(iAlpha, 2 * iAlpha + 1):
			term_U = 0.0j
			term_dUdT = 0.0j
			term_dUdxi = 0.0j
			term_d2Udxi2 = 0.0j
			term_d2UdT2 = 0.0j
			term_d2UdxidT = 0.0j
			aReal = float(iAlpha)
			for iJay in range(0, nn + iDelta + 1):
				jReal = float(iJay)
				gamma = mmReal / 2.0 + sReal + jReal

				term_U += (
					cADs[iAlpha][iDelta][iJay]
					* xi**iJay
					* bigT**(-gamma - 1.0 + iAlpha)
				)

				# These "more efficient" definitions seem to not work; numerical division issues?
				# term_dUdT += term_U*(-gamma-1.0+aReal)/bigT
				# term_d2UdT2 += term_dUdT*(-gamma-2.0+aReal)/bigT
				# term_dUdxi += term_U*jReal/xi
				# term_d2Udxi2 += term_dUdxi*(jReal-1.0)/xi
				# term_d2UdxidT += term_dUdT*jReal/xi

				# These explicit definitions work better.
				term_dUdT += (
					cADs[iAlpha][iDelta][iJay]
					* xi**iJay
					* bigT**(-gamma - 2.0 + aReal)
					* (-gamma - 1.0 + aReal)
				)
				term_d2UdT2 += (
					cADs[iAlpha][iDelta][iJay]
					* xi**iJay
					* bigT**(-gamma - 3.0 + aReal)
					* (-gamma - 1.0 + aReal)
					* (-gamma - 2.0 + aReal)
				)
				term_dUdxi += (
					cADs[iAlpha][iDelta][iJay]
					* jReal
					* xi**(iJay - 1)
					* bigT**(-gamma - 1.0 + aReal)
				)
				term_d2Udxi2 += (
					cADs[iAlpha][iDelta][iJay]
					* jReal
					* (jReal - 1.0)
					* xi**(iJay - 2)
					* bigT**(-gamma - 1.0 + aReal)
				)
				term_d2UdxidT += (
					cADs[iAlpha][iDelta][iJay]
					* jReal
					* xi**(iJay - 1)
					* bigT**(-gamma - 2.0 + aReal)
					* (-gamma - 1.0 + aReal)
				)
			# for i Jay
		ec2ba = (epsilonc2 / beta)**iAlpha
		ceb[0] += (term_U * ec2ba)
		dUdT += (term_dUdT * ec2ba)
		d2UdT2 += (term_d2UdT2 * ec2ba)
		dUdxi += (term_dUdxi * ec2ba)
		d2Udxi2 += (term_d2Udxi2 * ec2ba)
		dUdb += (term_U * (-iAlpha * ec2ba / beta))
		d2UdxidT += (term_d2UdxidT * ec2ba)
		d2Udxidb += (term_dUdxi * (-iAlpha * ec2ba / beta))
		d2UdbdT += (term_dUdT * (-iAlpha * ec2ba / beta))
		# for iDelta
	# for iAlpha

	ceb[0] *= lambda_nm
	dUdT *= lambda_nm
	d2UdT2 *= lambda_nm
	dUdxi *= lambda_nm
	d2Udxi2 *= lambda_nm
	dUdb *= lambda_nm
	d2UdxidT *= lambda_nm
	d2Udxidb *= lambda_nm
	d2UdbdT *= lambda_nm

	# calculating the time-domain EM fields
	bzr = beta * zr

	if polar == 2:  # radial polariation
		# E_rho, c
		ceb[1] = -(ii / rho) * (
			mmReal
			* ((nnReal + mmReal + 1.0) / bzr)
			* ceb[0]
			- (
				(2.0 * omega0 * xi / (sReal * zr))
				* d2UdbdT
			)
			+ (
				(omega0 / sReal)
				* (2.0 * xi * (nnReal + mmReal + 2.0) / bzr + mmReal * (xi / bzr + 1.0 / c))
				* dUdT
			)
			+ (xi * (2.0 * nnReal + 3.0 * mmReal + 4.0) / bzr) * dUdxi
			- (mmReal / zr) * dUdb
			+ (2.0 * xi**2 / bzr) * d2Udxi2
			+ (2.0 * omega0 * xi / sReal) * (2.0 * xi / bzr + 1.0 / c) * d2UdxidT
			- (2.0 * xi / zr) * d2Udxidb
			+ (2.0 * omega0**2 * xi / sReal**2) * (xi / bzr + 1.0 / c) * d2UdT2
		)

		# E_phi, c
		ceb[2] = (mmReal / rho) * (
			((nnReal + mmReal + 1.0) / bzr)
			* ceb[0]
			+ (omega0 / sReal) * (xi / bzr + 1.0 / c) * dUdT
			+ (xi / bzr) * dUdxi
			- (1.0 / zr) * dUdb
		)

		# E_z, c
		ceb[3] = (xi / rho**2) * (
			(-4.0 * omega0 / sReal) * (mmReal + 1.0) * dUdT
			- 4.0 * (mmReal + 1.0) * dUdxi
			- (4.0 * omega0**2 * xi / sReal**2) * d2UdT2
			- 4.0 * xi * d2Udxi2
			- (8.0 * omega0 * xi / sReal) * d2UdxidT
		)

		# B_rho, c
		ceb[4] = -(mmReal * omega0 / (c2 * sReal * rho)) * dUdT

		# B_phi
		ceb[5] = (
			-(ii * omega0 / (c2 * sReal * rho))
			* (mmReal * dUdT + 2.0 * xi * ((omega0 / sReal) * d2UdT2 + d2UdxidT))
		)

		# B_z = 0 for RP fields
		# ceb[6]=0.0j # already set as the default value

		return ceb
	else:
		print("Polarization not defined")
		exit()
	return
