import numpy as np
import sys 
from scipy import interpolate,optimize,integrate
import collections
import matplotlib.pyplot as plt


def OD2Biomass(od):
	return 0.0160348 + 0.218138*od

def convertCSVLineToFloats(line,delimiter):
	return [float(val) for val in line.rstrip().split(delimiter)]

def interpolateStandards(standardsFile,metabolites):
	interpolationFunctions = {}
	for meta in metabolites:
		matrix =[convertCSVLineToFloats(line,",") for line 
			in open(standardsFile+meta+".csv","rU").readlines()]
		matrix = [[0,0]]+[[row[0],np.mean(row[1:])] for row in matrix]
		if any(x<0 for x in [row[1] for row in matrix]): raise Exception("Bad Standard File")
		interpolationFunctions[meta] = interpolate.interp1d([row[1] for row in matrix],
			[row[0] for row in matrix],fill_value = "extrapolate")
	return interpolationFunctions

def computeTNeg1Flux(flowRate,standardsInterpolation,hplcValue,biomass,mw):
	return max([0,standardsInterpolation(hplcValue)*flowRate/biomass/mw])

def interpolateVolume(volumeFile,delimiter):
	matrix = [[float(x.rstrip().split(delimiter)[0]),float(x.rstrip().split(delimiter)[-1])]
		for x in open(volumeFile,"rU").readlines()[1:]]
	return interpolate.interp1d([row[0] for row in matrix],
			[row[1] for row in matrix],fill_value = "extrapolate")

def exponentialGrowthCurve(t,initialBMConc,c1):
	if isinstance(t,collections.Iterable):
		return [initialBMConc*(2**(c1*t1)) for t1 in t]
	else:
		return initialBMConc*(2**(c1*t))

def fitExpGrowthModelFromBiomassConc(odMatrix,plot=False):
	time = [row[0] for row in odMatrix]
	medianBiomassConcs = [OD2Biomass(np.median([np.mean(row[1:4]),np.mean(row[4:7]),np.mean(row[7:])])) for row in odMatrix]
	params, fit = optimize.curve_fit(exponentialGrowthCurve,time,medianBiomassConcs,
		bounds = (0,[max(medianBiomassConcs),1.0/60]))
	meanSError = mse(medianBiomassConcs,exponentialGrowthCurve(time,params[0],params[1]))
	if plot:
		plt.figure()
		plt.scatter(time,medianBiomassConcs)
		plt.xlabel("time")
		plt.ylabel("Conc Biomass")
		plt.title(str(meanSError))
		plt.plot(range(int(max(time))),exponentialGrowthCurve(range(int(max(time))),params[0],params[1]))

	return params,meanSError
def exponentialGrowthRate(initialBMConc,c1,t):
	return (2**(c1*t))*initialBMConc*c1*np.log(2)

def glucoseConc(t,c0,cm,t0,initBiomass,c1):
	glucoseConcIntegral = lambda t1: integrate.quad(glucoseConcIntegrand,0,t1,args=(c0,cm,t0,initBiomass,c1,))[0] 
	if isinstance(t,collections.Iterable):
		return [glucoseConcIntegral(t2) for t2 in t]
	else:
		return glucoseConcIntegral(t)

def glucoseConcIntegrand(t,c0,cm,t0,initBiomass,c1):
	return (c0*exponentialGrowthCurve(t,initBiomass,c1))/(1+2**(cm*(t0-t)))

def fitGlucoseModel(c1,initBiomass,glucoseConcMatrix,interpFunc,plot=False):
	time = [row[0] for row in glucoseConcMatrix]
	T0Conc = interpFunc(np.mean(glucoseConcMatrix[0][1:]))

	conc = [T0Conc - interpFunc(np.mean(row[1:])) for row in glucoseConcMatrix]
	partiallyParamerterizedFunc = lambda t,c0,cm,t0:glucoseConc(t,c0,cm,t0,initBiomass,c1)
	
	params,fit = optimize.curve_fit(partiallyParamerterizedFunc,time,conc,bounds = ([0,0,-1000],[1,1000,1000]))
	
	meanSError = mse(conc,glucoseConc(time,params[0],params[1],params[2],initBiomass,c1))
	if plot:
		plt.figure()
		plt.scatter(time,conc)
		plt.xlabel("time")
		plt.ylabel("Conc Glucose")
		plt.title(str(meanSError))
		plt.plot(range(int(max(time))),glucoseConc(range(int(max(time))),params[0],params[1],params[2],initBiomass,c1))

	return params,meanSError

def glucoseFlux(t,c0,cm,t0,initBiomass,c1,scale,mw):
	return gramsGlucoseConsumptionRatePerUnitBiomass(t,c0,cm,t0,initBiomass,c1)*scale/mw

def mse(measured,predicted):
	return np.mean(np.power(np.subtract(measured,predicted),2))

def gramsGlucoseConsumptionRatePerUnitBiomass(t,c0,cm,t0,initBiomass,c1):
	return c0/(1+2**(cm*(t-t0)))