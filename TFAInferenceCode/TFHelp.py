import sys
from gurobipy import *
import numpy as np
import pickle
import copy
from multiprocessing import Process, Queue
import time
import random

def TFVAMP(CS, GE, TFAOriginal, columns, q, iterLimit=20):
	"""
		Perform the variable TFA analysis on entire TFA matrix with complete GE matrix. The results for the
		range of possible values as well as the difference in factorization qualtiy are written to pickle files
		with the mean TFA activity written to a csv
	"""
	TFAVariable = list(TFAOriginal)
	TFVARange = [[[col, col] for col in row] for row in TFAOriginal]
	TFVAQuality = [[[-1, -1] for col in row] for row in TFVARange]
	for col in range(len(TFAOriginal[0])):
		TFA = [[row[col]] for row in TFAOriginal]
		singelGE = [[row[col]] for row in GE]
		tfRange, meanTFActivities, algorithmResult = TFAHelper.TFVASingleCondition(CS, singelGE, TFA, iterLimit)

		tempTFA = np.transpose(TFAVariable).tolist()
		tempTFA[col] = meanTFActivities
		TFAVariable = np.transpose(tempTFA).tolist()

		tempTFVARange = np.transpose(TFVARange, (1, 0, 2)).tolist()
		tempTFVARange[col] = tfRange
		TFVARange = np.transpose(tempTFVARange, (1, 0, 2)).tolist()

		tempTFVAQuality = np.transpose(TFVAQuality, (1, 0, 2)).tolist()
		tempTFVAQuality[col] = algorithmResult
		TFVAQuality = np.transpose(tempTFVAQuality, (1, 0, 2)).tolist()

	q.put([columns, TFAVariable, TFVARange, TFVAQuality])

	return 0


def TFVAMPRun(CS, GE, TFAOriginal, iterLimit=20, write=True, outfile="TFVAResult", numCores=2):
	colsPerProc = int(len(GE[0]) / numCores)
	q = Queue()
	startCol = 0
	processes = []
	TFAVariable = [[col for col in row] for row in TFAOriginal]
	TFARange = [[[-1, -1] for col in row] for row in TFAOriginal]
	TFAQuality = [[[-1, -1] for col in row] for row in TFAOriginal]

	for p in range(numCores):
		if p == numCores - 1:
			tempGE = [row[startCol:] for row in GE]
			tempTFA = [row[startCol:] for row in TFAOriginal]
			endCol = -1
		else:
			endCol = startCol + colsPerProc
			tempGE = [row[startCol:endCol] for row in GE]
			tempTFA = [row[startCol:endCol] for row in TFAOriginal]
		p = Process(target=TFVAMP, args=(CS, tempGE, tempTFA, [startCol, endCol], q))
		p.start()
		processes.append(p)
		startCol = endCol

	while True:
		running = any(p.is_alive() for p in processes)
		while not q.empty():
			[columns, TFAVar, TFVARan, TFVAQual] = q.get()
			if columns[1] == -1:
				TFAVariable = np.transpose(TFAVariable)
				TFAVariable[columns[0]:] = np.transpose(TFAVar)
				TFAVariable = np.transpose(TFAVariable).tolist()

				TFARange = np.transpose(TFARange, (1, 0, 2)).tolist()
				TFARange[columns[0]:] = np.transpose(TFVARan, (1, 0, 2))
				TFARange = np.transpose(TFARange, (1, 0, 2)).tolist()

				TFAQuality = np.transpose(TFAQuality, (1, 0, 2)).tolist()
				TFAQuality[columns[0]:] = np.transpose(TFVAQual, (1, 0, 2))
				TFAQuality = np.transpose(TFAQuality, (1, 0, 2)).tolist()
			else:
				TFAVariable = np.transpose(TFAVariable)
				TFAVariable[columns[0]:columns[1]] = np.transpose(TFAVar)
				TFAVariable = np.transpose(TFAVariable).tolist()

				TFARange = np.transpose(TFARange, (1, 0, 2)).tolist()
				TFARange[columns[0]:columns[1]] = np.transpose(TFVARan, (1, 0, 2))
				TFARange = np.transpose(TFARange, (1, 0, 2)).tolist()

				TFAQuality = np.transpose(TFAQuality, (1, 0, 2)).tolist()
				TFAQuality[columns[0]:columns[1]] = np.transpose(TFVAQual, (1, 0, 2))
				TFAQuality = np.transpose(TFAQuality, (1, 0, 2)).tolist()
		if not running:
			break
		

	if write:
		file = open(outfile + ".csv", "w")
		TFAHelper.writeMatrix2csv(TFAVariable, file)
		file.close()
		file = open(outfile + "Range.pkl", "wb")
		pickle.dump(TFARange, file)
		file.close()
		file =  open(outfile + "Quality.pkl", "wb")
		pickle.dump(TFAQuality,file)
		file.close()

	return TFAVariable, TFARange, TFAQuality


class TFAHelper():

	def __init__(self,TFAWholeMatrix, *args, **kwargs):
		self.common2SystematicName = {str.lower(x.rstrip().split(",")[0]):str.lower(x.rstrip().split(",")[1]) for x in open("YeastCommonAndSystematicGeneNames.csv","r").readlines()[1:]}
		self.systematic2Common = {value:key for key,value in self.common2SystematicName.items()}
		self.TFAWholeMatrix = copy.deepcopy(TFAWholeMatrix)
		pass

	@staticmethod
	def dense_optimize(c, A,b,rhsUB, lb, ub, vtype,
	                   solution,direction):
		
		#THIS CODE IS BASED OFF OF http://www.gurobi.com/documentation/8.0/examples/mip1_py.html
	  """
	  	Perform a maximization or minimization (based on the value of direction) of Ax (sense) rhs. Where sense determines
	  	the type of relation, less than/greater than/equal/etc. cx is the objective function and lb <= x <= ub. vtype determines
	  	the type of variables to use and solution holds the result
	  """
	  rows = len(A)
	  cols = len(A[0])

	  model = Model()

	  # Add variables to model
	  vars = []
	  for j in range(cols):
	    vars.append(model.addVar(lb=lb[j], ub=ub[j], vtype=vtype[j]))

	  # Populate A matrix

	  senseUB = [GRB.LESS_EQUAL for _ in rhsUB]


	  for i in range(rows):
	    expr = LinExpr()
	    for j in range(cols):
	      if A[i][j] != 0:
	        expr += A[i][j]*vars[j]
	    expr -= b[i][0]
	    model.addConstr(QuadExpr(expr*expr), senseUB[i], rhsUB[i])


	  # Populate objective
	  obj = LinExpr()

	  for j in range(cols):
	    if c[j] != 0:
	      obj += c[j]*vars[j]
	  model.setObjective(obj,direction)

	  # Solve
	  model.setParam("OutputFlag",False)
	  model.optimize()

	  # Write model to a file
	  #model.write('dense.lp')

	  if model.status == GRB.Status.OPTIMAL:
	    x = model.getAttr('x', vars)
	    for i in range(cols):
	      solution[i] = x[i]
	    return True
	  else:
	    return False

	@staticmethod
	def TFVASingleCondition(CS, GE, TFAOriginal, iterLimit=20):
		"""
			Perform TF activity variability analysis by sequentially maximimizing and minimizing
			each TF's activity individually subjust to the contraint that the TF activity solves
			the system with qualityCoefficient
		"""

		pGE = np.dot(CS, TFAOriginal).tolist()  # determine GE vector current facotrization solves
		meanError = TFAHelper.standardError(pGE, GE)  # compute error between current factorization and true expression

		lb = [0 for _ in range(len(CS[0]))]  # form lower bound for TF activity
		lb[-1] = TFAOriginal[-1][0] - TFAOriginal[-1][0] / 100  # fix pseudo TF activity
		ub = [GRB.INFINITY for _ in range(len(CS[0]))]  # network is not bounded above
		ub[-1] = TFAOriginal[-1][0] + TFAOriginal[-1][0] / 100  # fix upper bound for pseudo TF

		vtype = [GRB.CONTINUOUS for _ in range(len(CS[0]))]
		c = [0 for _ in range(len(CS[0]))]  # vector to hold which TF is objective
		tfRange = [[] for x in range(len(CS[0]))]  # will hold the range of values for the TFAs
		algorithmResult = [[] for x in
						   range(len(CS[0]))]  # will hold the relaxation of the current factorization result

		for tf in range(len(CS[0])):
			c[tf] = 1
			for direction in [GRB.MINIMIZE, GRB.MAXIMIZE]:
				qualityCoefficient = 0
				success = False
				iter = 0
				while not success and iter < iterLimit:
					qualityCoefficient += .1
					sol = [-1 for _ in range(len(CS[0]))]
					rhsUB = [qualityCoefficient for _ in range(len(CS))]
					success = TFAHelper.dense_optimize(c, CS, pGE, rhsUB, lb, ub, vtype, sol, direction)
					iter += 1
				if success:
					tfRange[tf].append(sol[tf])
					temp = np.dot(CS, TFAHelper.rowv2colv(sol)).tolist()
					tempError = TFAHelper.standardError(temp, GE)
					algorithmResult[tf].append(tempError - meanError)
				else:
					tfRange[tf].append(TFAOriginal[tf][0])
					algorithmResult[tf].append(-1)
			c[tf] = 0

		meanTFActivities = [np.mean(tf) for tf in tfRange]
		return tfRange, meanTFActivities, algorithmResult
	@staticmethod
	def TFVA(CS, GE, TFAOriginal, iterLimit=20, write=True, outfile="TFVAResult"):
		"""
			Perform the variable TFA analysis on entire TFA matrix with complete GE matrix. The results for the
			range of possible values as well as the difference in factorization qualtiy are written to pickle files
			with the mean TFA activity written to a csv
		"""

		TFAVariable = list(TFAOriginal)
		TFVARange = [[[col, col] for col in row] for row in TFAOriginal]
		TFVAQuality = [[[-1, -1] for col in row] for row in TFVARange]
		for col in range(len(TFAOriginal[0])):
			TFA = [[row[col]] for row in TFAOriginal]
			singelGE = [[row[col]] for row in GE]
			tfRange, meanTFActivities, algorithmResult = TFAHelper.TFVASingleCondition(CS, singelGE, TFA, iterLimit)

			tempTFA = np.transpose(TFAVariable).tolist()
			tempTFA[col] = meanTFActivities
			TFAVariable = np.transpose(tempTFA).tolist()

			tempTFVARange = np.transpose(TFVARange, (1, 0, 2)).tolist()
			tempTFVARange[col] = tfRange
			TFVARange = np.transpose(tempTFVARange, (1, 0, 2)).tolist()

			tempTFVAQuality = np.transpose(TFVAQuality, (1, 0, 2)).tolist()
			tempTFVAQuality[col] = algorithmResult
			TFVAQuality = np.transpose(tempTFVAQuality, (1, 0, 2)).tolist()

		if write:
			TFAHelper.writeMatrix2csv(TFAVariable, open(outfile + ".csv", "w"))
			pickle.dump(TFVARange, open(outfile + "Range.pkl", "wb"))
			pickle.dump(TFVAQuality, open(outfile + "Quality.pkl", "wb"))

		return TFAVariable, TFVARange, TFVAQuality

	@staticmethod
	def rowv2colv(a):
			"""
			Transpose a 1D vector
			"""
			try: 
				a[0][0]
				return [e[0] for e in a]
			except:
				return [[e] for e in a]

	@staticmethod
	def standardError(a,b):
		"""
		Compute standard Error
		"""
		return np.mean([abs(x[0]-y[0]) for x,y in zip(a,b)])

	@staticmethod
	def writeMatrix2csv(matrix,file):
		"""
		Write a nested list of list to a csv file
		"""
		file.write(str(matrix[0][0]))
		[file.write(","+str(val)) for val in matrix[0][1:]]
		for row in matrix[1:]:
			file.write("\n"+str(row[0]))
			for val in row[1:]:
				file.write(","+str(val))

	def standardScore(self,matrix):
		return [[(matrix[row][col]-np.mean(self.TFAWholeMatrix[row]))/np.std(self.TFAWholeMatrix[row]) for col in range(len(matrix[row]))] for row in range(len(matrix))]

	def absoluteLog2FoldChange(self,matrix):
		return [[np.abs(np.log2(matrix[row][col]/matrix[row][-1])) for col in range(len(matrix[row]))] for row in range(len(matrix))]

	def rankingChange(self,matrix):
		tempMat = np.transpose(self.absoluteLog2FoldChange(matrix))
		for x in range(len(tempMat)):
			temp = list(tempMat[x])
			indices = range(len(temp))
			indices.sort(key=lambda x:temp[x],reverse=True)
			tempMat[x] = [indices.index(tf) for tf in range(len(temp))]
		return np.transpose(tempMat).tolist()


	def rankingChangePercentile(self,matrix):
		tempMat = np.transpose(self.absoluteLog2FoldChange(matrix))
		for x in range(len(tempMat)):
			tempMin = 0#min(tempMat[x])
			tempMax = max(tempMat[x]) - tempMin
			tempMat[x] = [(val-tempMin)/tempMax for val in tempMat[x]]
		return np.transpose(tempMat).tolist()

	def valueRelative2Max(self,matrix):
		tempMat = np.transpose(matrix)
		for x in range(len(tempMat)):
			temp = list(tempMat[x])
			temp = [(temp[i]/matrix[i][-1]) for i in range(len(temp))]
			tempMat[x] = [i/max(temp) for i in temp]
			return np.transpose(tempMat).tolist()

	def getPValuesForInteractions(self,interactions,allScores):
		allScores.sort(reverse=True)
		number = len(allScores)
		thresh = 1e-6
		interactionPvals = dict(interactions)
		for inter in interactions:
			index = 0.0
			for score in allScores:
				if interactions[inter] > score:
					index -= 1
					break
				index += 1
			interactionPvals[inter] = index/number
		return interactionPvals

	@staticmethod
	def flatten(l):
		return [item for sublist in l for item in sublist]

	def catagorizeInteractionsForSignificance(self,alpha,interactions,tailed=True):
		if tailed:
			significantFoundInteractions = [inter for inter in interactions if interactions[inter] < alpha/2 or interactions[inter] > (1-alpha/2)]
		else:
			significantFoundInteractions = [inter for inter in interactions if interactions[inter] < alpha]
		insignificantInteractions= [inter for inter in interactions if inter not in significantFoundInteractions]
		return significantFoundInteractions,insignificantInteractions

	def generateRandomMatrix(self,matrix):
		allScores = self.flatten(matrix)
		random.shuffle(allScores)
		colLength = len(matrix[0])
		rowLength = len(matrix)
		return [[allScores[row*colLength+col] for col in range(colLength)] for row in range(rowLength)]
		#return np.random.normal(np.mean(self.flatten(self.TFAWholeMatrix)),np.std(self.flatten(self.TFAWholeMatrix)),(len(matrix),len(matrix[0]))).tolist()


	def expectedInteractionFallInRankingLevel(self,interactionsFound,numberOfTFs):
		xcoords = []
		ycoords = []
		for tf in range(numberOfTFs):
			count = 0
			for inter in interactionsFound:
				if interactionsFound[inter] <= tf:
					count += 1
			xcoords.append(float(tf)/numberOfTFs)
			ycoords.append(float(count)/len(interactionsFound))
		return xcoords,ycoords

	def expectedInteractionFallInPercentileLevel(self,interactionsFound):
		xcoords = []
		ycoords = []
		for score in [x/100.0 for x in range(100)]:
			count = 0
			for inter in interactionsFound:
				if interactionsFound[inter] >= score:
					count += 1
			xcoords.append(score)
			ycoords.append(float(count)/len(interactionsFound))
		return xcoords,ycoords

	@staticmethod
	def dist(A,B,tol=5):
		A = [np.round(a,tol) for a in A]
		B = [np.round(b,tol) for b in B]

		def sgn(a):
			if a > 0: return 1
			elif a < 0: return -1
			else: return 0

		mA = np.round(np.mean(A),tol)
		mB = np.round(np.mean(B),tol)
		rA = max(A) - mA
		rB = max(B) - mB

		if rB+rA == 0:
			return mB - mA + 2 * sgn(mB-mA)
		elif mA != mB:
			return (mB - mA)/(rB+rA) + sgn(mB-mA)
		else:
			return (rB - rA)/max([rB,rA])
		# A = [np.round(a,tol) for a in A]
		# B = [np.round(b,tol) for b in B]

		# return min([abs(np.log2(max(A)/max(B))),abs(np.log2(min(A)/min(B)))])

		
		def sgn(a):
			if a > 0: return 1
			elif a < 0: return -1
			else: return 0

		mA = np.round(np.mean(A),tol)
		mB = np.round(np.mean(B),tol)
		rA = max(A) - mA
		rB = max(B) - mB

		if rB+rA == 0:
			return mB - mA + 2 * sgn(mB-mA)
		elif mA != mB:
			return (mB - mA)/(rB+rA) + sgn(mB-mA)
		else:
			return (rB - rA)/max([rB,rA])


	def rankingChangeV(self,varMatrix):
		tempMat = np.transpose(varMatrix,(1,0,2)).tolist()
		for x in range(len(tempMat)):
			temp = list(tempMat[x])
			temp = [abs(TFAHelper.dist(temp[i],varMatrix[i][-1])) for i in range(len(temp))]
			indices = range(len(temp))
			indices.sort(key=lambda x:temp[x],reverse=True)
			tempMat[x] = [indices.index(tf) for tf in range(len(temp))]
		return np.transpose(tempMat).tolist()

	@staticmethod
	def distFromWTRange(varMatrix):
		tempMat = copy.deepcopy(varMatrix)
		for x in range(len(varMatrix)):
			tempMat[x] = [np.abs(TFAHelper.dist(tempMat[x][col],tempMat[x][-1])) for col in range(len(varMatrix[x]))]
		return tempMat

	@staticmethod
	def transformScoredMatrixToRank(matrix):
		tempMat = np.transpose(matrix).tolist()
		for x in range(len(tempMat)):
			temp = list(tempMat[x])
			indices = range(len(temp))
			indices.sort(key=lambda x:temp[x],reverse=True)
			tempMat[x] = indices
		return np.transpose(tempMat).tolist()

	@staticmethod
	def BinDistMatrixV(distMatrix):
		def discrete(val):
			if val < .01:
				return 0
			elif val <= 1:
				return 1
			elif val <= 2:
				return 2
			else :
				return 3
		return [[discrete(col) for col in row] for row in distMatrix]
