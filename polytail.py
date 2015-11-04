#!/usr/bin/env python
# -*- coding: utf-8 -*-

class Tail:
	def __init__(self, ylist, i, dh, x0, xc):
		self.x0 = x0
		self.xc = xc

		self.f   = ylist[i-1]
		self.df  = self.deriv(ylist, i-1, dh, 1)
		self.ddf = self.deriv(ylist, i-1, dh, 2)

		self.colist = self.coefficients(self.f, self.df, self.ddf, self.x0, self.xc)

	def deriv(self, ylist, i0, dh, order):
		if order == 1:
			return  (1/12.*ylist[i0-2] - 2/3.*ylist[i0-1] + 2/3.*ylist[i0+1] - 1/12.*ylist[i0+2])/dh

		if order == 2:
			return (-1/12.*ylist[i0-2] + 4/3.*ylist[i0-1] -5/2.*ylist[i0] + 4/3.*ylist[i0+1] - 1/12.*ylist[i0+2])/dh/dh

	def coefficients(self, f, df, ddf, x0, xc):
		colist = [0,0,0,0,0,0]

		colist[0] = -(xc**3 * (2*f * (10*x0**2 - 5*x0*xc + xc**2) + x0 * (x0-xc) * (ddf*x0 * (x0 - xc) + 2*df * (-4*x0+xc)))) / (2.*(x0 - xc)**5)
		colist[1] = (xc**2 * (60*f * x0**2 - (x0 - xc) * (2*df * (6*x0 - xc) * (2*x0 + xc) + ddf*x0 * (-3*x0**2 + x0*xc + 2*xc**2)))) / (2.*(x0-xc)**5)
		colist[2] = (xc * (-60*f * x0*(x0 + xc) + (x0 - xc) * (12*df * x0*(2*x0 + 3*xc) - ddf*(x0 - xc) * (3*x0**2 + 6*x0*xc+xc**2)))) / (2.*(x0-xc)**5)
		colist[3] = (20*f * (x0**2 + 4*x0*xc + xc**2) + (x0 - xc) * (ddf * (x0 - xc) * (x0**2 + 6*x0*xc + 3*xc**2) - 4*df*(2*x0**2 + 10*x0*xc + 3*xc**2))) / (2.*(x0-xc)**5)
		colist[4] = (-30*f * (x0 + xc) + (x0 - xc) * (-ddf * (x0 - xc) * (2*x0 + 3*xc) + 2*df*(7*x0 + 8*xc))) / (2.*(x0-xc)**5)
		colist[5] = (12*f + (-6*df + ddf*(x0 - xc)) * (x0 - xc)) / (2.*(x0-xc)**5)

		return colist

	def tail(self, x):
		if (x < self.x0) or (x > self.xc):
			return 0.0

		return self.colist[0] + self.colist[1]*x + self.colist[2]*x**2 + self.colist[3]*x**3 + self.colist[4]*x**4 + self.colist[5]*x**5
