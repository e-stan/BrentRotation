import numpy as np
import sys 
from scipy import interpolate


def OD2Biomass(od):
	return 0.0160348 + 0.218138*od

def convertCSVLineToFloats(line,delimiter):
	return [float(val) for val in line.split(delimiter)]

def interpolateStandards(standardsFile,metabolites):
	interpolationFunctions = {}
	for meta in metabolites:
		matrix =[convertCSVLineToFloats(line,",") for line 
			in open(standardsFile+meta+".csv","r").readlines()]
		matrix = [[0,0]]+[[row[0],np.mean(row[1:])] for row in matrix]	
		interpolationFunctions[meta] = interpolate.interp1d([row[1] for row in matrix],
			[row[0] for row in matrix],fill_value = "extrapolate")
	return interpolationFunctions

def computeTNeg1Flux(flowRate,standardsInterpolation,hplcValue,biomass,mw):
	return max([0,standardsInterpolation(hplcValue)*flowRate/biomass/mw])

def interpolateVolume(volumeFile,delimiter):
	matrix = [[float(x.rstrip().split(delimiter)[0]),float(x.rstrip().split(delimiter)[-1])]
		for x in open(volumeFile,"r").readlines()[1:]]
	return interpolate.interp1d([row[0] for row in matrix],
			[row[1] for row in matrix],fill_value = "extrapolate")

def exponentialGrowthCurve(t,initialBMConc,c1):
	return initialBMConc*(2**(c1*time))

def fitExpGrowthModelFromBiomassConc(odMatrix):
	time = [row[0] for row in odMatrixd]
	medianODs = np.me
	return scipy.optimize.curve_fit(exponentialGrowthCurve,time,biomass,
		bounds = (0,[biomass[-1],1/60]))

def exponentialGrowthRate(initialBMConc,c1,t):
	return (2**(c1*t))*initialBMConc*c1*np.log(2)