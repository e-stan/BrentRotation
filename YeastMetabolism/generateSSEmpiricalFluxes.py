import sys 
from chemostatHelper import *

directory = sys.argv[1]
outputFile = sys.argv[2]
strain = sys.argv[3]
accuracy = 5

relevantMetabolites = ["aceticAcid","ethanol","glycerol"]

molarMasses = [60.05,46.068,92.0938]

writeHeaders = False

try:
	checkOutput = open(outputFile,"r")
	checkOutput.close()
	outputFile = open(outputFile,"a")
except:
	outputFile = open(outputFile,"w")
	outputFile.write("strain,")
	for meta in relevantMetabolites:
		outputFile.write(meta+",")
	outputFile.write("growthRate")

hplcStandards = "hplcStandards/"


metaboliteMolecularWeights = {meta:val for meta,val in zip(relevantMetabolites,
	molarMasses)}

odFile = "od.csv"

try:
	interpFunction = interpolateStandards(directory+hplcStandards,relevantMetabolites)

	timeScale = 60
	flowRate = float(open(directory+"mediaFlowRate.csv","r").readline())*timeScale

	odFile = open(directory+odFile,"r")

	# od = convertCSVLineToFloats(odFile.readline(),",")[1:]
	# od = np.median([np.mean(od[0:3]),np.mean(od[3:6]),np.mean(od[6:])])

	volumeInterpolator = interpolateVolume(directory+"volume.csv",",")

	#biomass = OD2Biomass(od)*volumeInterpolator(-1)

	biomassVolume = .03
	biomassConc = np.mean(convertCSVLineToFloats(open(directory+"dryWeight.csv","r").readline(),",")[1:])/biomassVolume
	biomass = biomassConc*volumeInterpolator(-1)

	metaFluxValues = []
	for meta in relevantMetabolites:
		metaFluxValues.append(np.round(np.mean(
			[computeTNeg1Flux(flowRate,interpFunction[meta],val,biomass,metaboliteMolecularWeights[meta]) for val
			in convertCSVLineToFloats(open(directory+meta+".csv","r").readline(),",")[1:]]),accuracy))

	odT0 = convertCSVLineToFloats(odFile.readline(),",")[1:]
	odT0 = np.median([np.mean(odT0[0:3]),np.mean(odT0[3:6]),np.mean(odT0[6:])])
	biomassT0 = OD2Biomass(odT0)

	timeT0ToTNeg1 = -1*float(open(directory+"timeTNeg1ToT0.csv","r").readline().rstrip().split(",")[0])
	growthRate = (biomassT0-biomassConc)*volumeInterpolator(-1)*(timeScale/timeT0ToTNeg1)/biomass
	
	outputFile.write("\n"+strain)
	for meta in metaFluxValues: outputFile.write(","+str(meta))
	outputFile.write(","+str(np.round(growthRate,accuracy)))
    
except OSError as e:
	pass
except:
	print("error")





