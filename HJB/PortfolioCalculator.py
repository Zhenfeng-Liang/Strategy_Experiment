from Model import Model
import numpy as np
from scipy import linalg

class PortfolioCalculator:

	def __init__(self, model, corrMatr):
		self.model = model
		self.corrMatr = corrMatr

	def instCov(self, x):
		
		dofx = self.model.diffV(x)
		covx = dofx.dot(self.corrMatr).dot(dofx.T)

		return covx

	def instCovDer(self, x, i):

		dofx = self.model.diffV(x)
		ddofx = self.model.diffDer(x, i)
		derx = ddofx.dot(self.corrMatr).dot(dofx.T) + dofx.dot(self.corrMatr).dot(ddofx.T)

		return derx

	def instCovDer2(self, x, i, j):

		dofx = self.model.diffV(x)
		ddofxi = self.model.diffDer(x, i)
		ddofxj = self.model.diffDer(x, j)
		d2dofxij = self.model.diffDer2(x, i, j)

		derx = d2dofxij.dot(self.corrMatr).dot(dofx.T) \
		+ ddofxi.dot(self.corrMatr).dot(ddofxj.T) \
		+ ddofxj.dot(self.corrMatr).dot(ddofxi.T) \
		+ dofx.dot(self.corrMatr).dot(d2dofxij.T)

		return derx

	def invInstCov(self, x):

		invCovx = linalg.pinv(self.instCov(x))

		return invCovx

	def invInstCovDer(self, x, i):

		invCovx = self.invInstCov(x)
		covxdi = self.instCovDer(x, i)
		derx = -1.0 * invCovx.dot(covxdi).dot(invCovx)

		return derx

	def invInstCovDer2(self, x, i, j):

		invCovx = self.invCovx(x)
		covxdi = self.instCovDer(x, i)
		covxdj = self.instCovDer(x, j)

		tmp = covxdi.dot(invCovx).dot(covxdj) + covxdj.dot(invCovx).dot(covxdi) - self.instCovDer2(x, i, j)
		derx = invCovx.dot(tmp).dot(invCovx)

		return derx