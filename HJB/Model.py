import numpy as np 

class Model:
    
    def __init__(self, param):
        self.modelType = param["modelType"]
        self.mu = param["mu"]
        self.vol = param["vol"]

        if "lambda" in param.keys():
            self.lam = param["lambda"]

    def driftV(self, x):

    	if self.modelType == "MeanReverting" or self.modelType == "CIR":
    		res = self.lam *  (self.mu - x)
    	elif self.modelType == "LogNormal":
    		res = self.mu * x

    	return res

    def driftDer(self, x, i):

		F = len(x)
		der = np.zeros(F)

    	if self.modelType == "MeanReverting" or self.modelType == "CIR":
    		der[i] = -self.lam[i]
    	elif self.modelType == "LogNormal":
			der[i] = self.mu[i]    		 
    	
    	return der

    def driftDer2(self, x, i, j):

    	F = len(x)
    	der2 = np.zeros(F)

    	return der2

    def diffV(self, x):

    	if self.modelType == "MeanReverting":
    		b = np.diag(self.vol)
    	elif self.modelType == "LogNormal":
    		b = np.diag(self.vol * x)
    	elif self.modelType == "CIR":
    		b = np.diag(self.vol * np.sqrt(x))

    	return b

    def diffDer(self, x, i):

    	
    	# Assuming only one random driver for each asset
    	F = len(x)
    	der = np.zeros((F, F))

    	if self.modelType == "MeanReverting":
    		# Do nothing
    	elif self.modelType == "LogNormal":
    		der[i][i] = self.vol[i]
    	elif self.modelType == "CIR":
    		der[i][i] = 0.5 * self.vol[i] * (x[i]**-0.5)

    	return der

    def diffDer2(self, x, i, j):

    	F = len(x)
    	der2 = np.zeros((F, F))

    	if self.modelType == "CIR":

    		if i == j:
    			der2[i][i] = -0.25 * self.vol[i] * (x[i] ** -1.5)
    	else:
    		# Do nothing for mean reverting and lognormal model

    	return der2
