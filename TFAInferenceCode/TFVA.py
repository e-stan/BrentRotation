import sys
from gurobipy import *
import numpy as np
import pickle

#THIS CODE IS BASED OFF OF http://www.gurobi.com/documentation/8.0/examples/mip1_py.html

def dense_optimize(c, A,b,rhsUB, lb, ub, vtype,
                   solution,direction):

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

def rowv2colv(a):
	try: 
		a[0][0]
		return [e[0] for e in a]
	except:
		return [[e] for e in a]

def standardError(a,b):
	return np.mean([abs(x[0]-y[0]) for x,y in zip(a,b)])

def writeMatrix2csv(matrix,file):
	file.write(str(matrix[0][0]))
	[file.write(","+str(val)) for val in matrix[0][1:]]
	for row in matrix[1:]:
		file.write("\n"+str(row[0]))
		for val in row[1:]:
			file.write(","+str(val))


def TFVASingleCondition(CS,GE,TFAOriginal,iterLimit = 20):
	"""
		Perform TF activity variability analysis by sequentially maximimizing and minimizing
		each TF's activity individually subjust to the contraint that the TF activity solves 
		the system with qualityCoefficient
	"""

	pGE = np.dot(CS,TFAOriginal).tolist()
	meanError = standardError(pGE,GE)

	lb = [0 for _ in range(len(CS[0]))]
	lb[-1] = TFAOriginal[-1][0]-TFAOriginal[-1][0]/100
	ub = [GRB.INFINITY for _ in range(len(CS[0]))]
	ub[-1] = TFAOriginal[-1][0]+TFAOriginal[-1][0]/100

	vtype = [GRB.CONTINUOUS for _ in range(len(CS[0]))]
	c = [0 for _ in range(len( CS[0]))]
	tfRange = [[] for x in range(len(CS[0]))]
	algorithmResult = [[] for x in range(len(CS[0]))]

	for tf in range(len(CS[0])):
		c[tf] = 1
		for direction in [GRB.MINIMIZE,GRB.MAXIMIZE]:
			qualityCoefficient = 0
			success = False
			iter = 0
			while not success and iter < iterLimit: 
				qualityCoefficient+=.1
				sol = [-1 for _ in range(len(CS[0]))]
				rhsUB = [qualityCoefficient for _ in range(len(CS))]
				success = dense_optimize(c,CS,pGE,rhsUB,lb,ub,vtype,sol,direction)
				iter+=1
			if success: 
				tfRange[tf].append(sol[tf])
				temp = np.dot(CS,rowv2colv(sol)).tolist()
				tempError = standardError(temp,GE)
				algorithmResult[tf].append(tempError-meanError)
			else:
				tfRange[tf].append(TFAOriginal[tf][0])
				algorithmResult[tf].append(-1)
		c[tf] = 0

	meanTFActivities = [np.mean(tf) for tf in tfRange]
	return tfRange,meanTFActivities,algorithmResult


def TFVA(CS,GE,TFAOriginal,iterLimit=20,write=True,outfile = "TFVAResult"):
	"""
		Perform the variable TFA analysis on entire TFA matrix with complete GE matrix. The results for the
		range of possible values as well as the difference in factorization qualtiy are written to pickle files
		with the mean TFA activity written to a csv
	"""

	TFAVariable = list(TFAOriginal)
	TFVARange = [[[col,col] for col in row] for row in TFAOriginal]
	TFVAQuality = [[[-1,-1] for col in row] for row in TFVARange]
	for col in range(len(TFAOriginal[0])):
		TFA = [[row[col]] for row in TFAOriginal]
		singelGE = [[row[col]] for row in GE]
		tfRange,meanTFActivities,algorithmResult = TFVASingleCondition(CS,singelGE,TFA,iterLimit)

		tempTFA = np.transpose(TFAVariable).tolist()
		tempTFA[col] = meanTFActivities
		TFAVariable = np.transpose(tempTFA).tolist()

		tempTFVARange = np.transpose(TFVARange,(1, 0, 2)).tolist()
		tempTFVARange[col] = tfRange
		TFVARange = np.transpose(tempTFVARange,(1, 0, 2)).tolist()

		tempTFVAQuality = np.transpose(TFVAQuality,(1, 0, 2)).tolist()
		tempTFVAQuality[col] = algorithmResult
		TFVAQuality = np.transpose(tempTFVAQuality,(1, 0, 2)).tolist()

	if write:
		writeMatrix2csv(TFAVariable,open(outfile+".csv","w"))
		pickle.dump(TFVARange,open(outfile+"Range.pkl","wb"))
		pickle.dump(TFVAQuality,open(outfile+"Quality.pkl","wb"))


	return TFAVariable,TFVARange,TFVAQuality




