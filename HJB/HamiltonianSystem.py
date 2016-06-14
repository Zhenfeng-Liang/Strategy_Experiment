import numpy as np
from scipy import linalg

class HamiltonianSystem:

	def __init__(self, portCalc, utiCalc):

		self.portCalc = portCalc
		self.utiCalc = utiCalc
		self.oneOverGamma = oneOverGamma
		self.kappa = utiCalc.kappa

	def ham(self, x, p):

		a = self.portCalc.model.driftV(x)

		term1 = 0.5 * p.T.dot(self.portCalc.instCov(x)).dot(p)
		term2 = self.oneOverGamma * p.T.dot(a)
		term3 = 0.5 * self.kappa * self.oneOverGamma * a.T.dot(\
			linalg.solve(self.portCalc.instCov(x), a))

		val = term1 + term2 + term3

		return val

	def hamDx(self, x, p):

		a = self.portCalc.model.driftV(x)
		covx = self.portCalc.instCov(x)

		F = len(x)
		val = np.zeros((F, 1))

		for i in range(F):

			term1 = 0.5 * p.T.dot(self.portCalc.instCovDer(x, i)).dot(p)
			term2 = self.oneOverGamma * p.T.dot(self.portCalc.model.driftDer(x,i))
			term3 = self.kappa * self.oneOverGamma * a.T.dot(\
				linalg.solve(covx, self.portCalc.model.driftDer(x, i)) + 0.5 * self.portCalc.invInstCovDer(x, i).dot(a))

			val[i][0] = term1 + term2 + term3

		return val


	def hamDxx(self, x, p):

		a = self.portCalc.model.driftV(x)
		covx = self.portCalc.instCov(x)

		F = len(x)
		val = np.zeros((F, F))

		for i in range(F):

			daxi = self.portCalc.model.driftDer(x, i)
			dinvCxi = self.portCalc.invInstCovDer(x, i)

			for j in range(F):

				daxj = self.portCalc.model.driftDer(x, j)
				d2axij = self.portCalc.model.driftDer2(x, i, j)
				dinvCxj = self.portCalc.invInstCovDer(x, j)
				d2invCxij = self.portCalc.invInstCovDer2(x, i, j)

				d2Cxij = self.portCalc.instCovDer2(x, i, j)

				term1 = 0.5 * p.T.dot(d2Cxij).dot(p)

				term2 = self.oneOverGamma * p.T.dot(d2axij)

				term3 = self.kappa * self.oneOverGamma * (\
					a.T.dot(linalg.solve(covx, d2axij)) +\
					daxj.T.dot(linalg.solve(covx, daxi)) +\
					a.T.dot(dinvCxj).dot(daxi) +\
					a.T.dot(dinvCxi).dot(daxj) + \
					0.5 * a.T.dot(d2invCxij).dot(a))

				val[i][j] = term1 + term2 + term3

		return val


	def hamDp(self, x, p):

		val = self.portCalc.instCov(x).dot(p) +\
		 self.oneOverGamma * self.portCalc.model.driftV(x)

		return val

	def hamDpp(self, x, p):

		val = self.portCalc.instCov(x)

		return val

	def hamDxp(self, x, p):

		F = len(x)
		val = np.zeros((F, F))

		for i in range(F):
			val[i][:] = (self.portCalc.instCovDer(x, i).dot(p) +\
			self.oneOverGamma * self.portCalc.model.driftDer(x, i)).T

		return val
