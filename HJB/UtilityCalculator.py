from math import exp

class UtilityCalculator:

	def __init__(self, gamma, uType, a = 1.0, b = 1.0):

		self.a = a
		self.b = b
		self.gamma = gamma
		self.oneOverGamma = 1.0 / gamma
		self.uType = uType

		if self.uType == "HARA" or self.uType == "CRRA":
			self.kappa = self.oneOverGamma - 1.0
		elif self.uType == "CARA":
			self.kappa = -1.0

	def Au(self, v):

		if self.uType == "CARA":
			res = self.gamma
		elif self.uType == "CRRA":
			res = self.gamma / v
		elif self.uType == "HARA":
			res = self.b / (self.a + self.b / self.gamma * v) 

		return res

	def U(self, v):

		if self.uType == "CARA":
			res = -exp(-self.gamma * v) / self.gamma
		elif self.uType == "CRRA":
			res = v ** (1.0 - self.gamma) / (1.0 - self.gamma)
		elif self.uType == "HARA":
			res = self.gamma / (1.0 - self.gamma) \
			* (self.a + self.b / self.gamma * v) ** (1.0 - self.gamma)

		return res

	def UDer(self, v):

		if self.uType == "CARA":
			res = exp(-self.gamma * v)
		elif self.uType == "CRRA":
			res = v ** (-self.gamma)
		elif self.uType == "HARA":
			res = self.b * (self.a + self.b / self.gamma * v) ** (-self.gamma)

		return res

	def UDer2(self, v):

		if self.uType == "CARA":
			res = -self.gamma * exp(-self.gamma * v)
		elif self.uType == "CRRA":
			res = -self.gamma * v ** (-self.gamma - 1.0)
		elif self.uType == "HARA":
			res = -self.b ** 2 \
			* (self.a + self.b / self.gamma * v) ** (-self.gamma - 1.0)

		return res
		