import sys 
from chemostatHelper import *
import matplotlib.pyplot as plt

print(glucoseConc(150,.03333,35.5151,-14.4463,.560,.00506))

directoryBase = sys.argv[1]
outputFile = sys.argv[2]

relevantMetabolites = ["aceticAcid","ethanol","glycerol","glucose"]
molarMasses = [60.05,46.068,92.0938,180.156]
accuracy = 5
timeScale = 60
biomassVolume = .03
timePointsOfInterestInBolus = [30,90,300]


hplcStandards = "hplcStandards/"
odFile = "od.csv"
volumeFile = "volume.csv"
mediaFlowRateFile = "mediaFlowRate.csv"
dryWeightFile = "dryWeight.csv"
timeT0ToTNeg1File = "timeTNeg1ToT0.csv"


metadata =  [x.rstrip().split(",") for x in open(directoryBase+"metadata.csv","rU").readlines()]
headers = metadata[0]
metadata = metadata[1:-1]

headerToIndex = {header:pos for header,pos in zip(headers,range(len(headers)))}

outputFile = open(outputFile,"w")
outputFile.write("strain,")
for meta in relevantMetabolites[:4]:
	outputFile.write("TNeg1Flux"+meta+",")
outputFile.write("biomassTNeg1"+",")
outputFile.write("growthRateTNeg1"+",")
outputFile.write("ExponentialGrowthModelMSE,")
outputFile.write("glucoseModelFit,")
[outputFile.write("glucoseFlux"+str(time)+",") for time in timePointsOfInterestInBolus]

	
metaboliteMolecularWeights = {meta:val for meta,val in zip(relevantMetabolites,
	molarMasses)}

badFileBeyondJustNE = 0
goodFiles = 0

for sample in metadata:
	directory = directoryBase+sample[-1]+"/"
	strain = sample[headerToIndex["strainID"]]+"_"+sample[headerToIndex["date"]]+"_"+sample[headerToIndex["chemostat"]]
	if sample[headerToIndex["goodOD"]] == "yes" and sample[headerToIndex["goodDW"]] == "yes" and sample[headerToIndex["goodHPLC"]] == "yes":
		goodFiles+=1
		try:
			if goodFiles < 10:
				plot = False
			else:
				plot = False

			interpFunction = interpolateStandards(directory+hplcStandards,relevantMetabolites)

			flowRate = float(open(directory+mediaFlowRateFile,"rU").readline())*timeScale

			ods = open(directory+odFile,"rU")
			ods.readline()

			volumeInterpolator = interpolateVolume(directory+volumeFile,",")
			biomassConc = np.median(convertCSVLineToFloats(open(directory+dryWeightFile,"rU").readline(),",")[1:])/biomassVolume
			biomass = biomassConc*volumeInterpolator(-1)

			metaFluxValues = []
			metaFiles = {}
			for meta in relevantMetabolites[:4]:
				metaFile = open(directory+meta+".csv","rU")
				metaFluxValues.append(np.round(np.mean(
					[computeTNeg1Flux(flowRate,interpFunction[meta],val,biomass,metaboliteMolecularWeights[meta]) for val
					in convertCSVLineToFloats(metaFile.readline(),",")[1:]]),accuracy))
				metaFiles[meta] = metaFile


			outputFile.write("\n"+strain)
			for meta in metaFluxValues: outputFile.write(","+str(meta))
			outputFile.write(","+str(biomass
				))
			params,fit = fitExpGrowthModelFromBiomassConc([convertCSVLineToFloats(line,",") for line in ods.readlines()],plot=plot)

			biomassT0 = params[0]

			c1 = params[1]

			timeT0ToTNeg1 = abs(float(open(directory+timeT0ToTNeg1File,"rU").readline().rstrip()))
			
			volume = volumeInterpolator(-1)

			#growthRate = (biomassT0-biomassConc)*(volume)*(timeScale/timeT0ToTNeg1)/biomass
			
			growthRate = exponentialGrowthRate(params[0],params[1],0)*timeScale/params[0]

			outputFile.write(","+str(np.round(growthRate,accuracy))+","+str(fit))

			#compute glucose model fit
	
			params,fit = fitGlucoseModel(params[1],params[0],[convertCSVLineToFloats(line,",") for line in open(directory+"glucose.csv","rU").readlines()[1:]],interpFunc = interpFunction["glucose"],plot=plot)

			outputFile.write(","+str(fit))

			[outputFile.write(","+str(flux)) for flux in [glucoseFlux(t,params[0],params[1],params[2],biomassT0,c1,timeScale,metaboliteMolecularWeights["glucose"]) for t in timePointsOfInterestInBolus]]




		
		except (OSError,IOError):
			#print(directory, "missing files")
			pass
		except Exception as inst:
			#print(directory, inst.args[0])
			badFileBeyondJustNE+=1
			pass
		except:
			#print(directory, sys.exc_info()[0])
			badFileBeyondJustNE+=1

print(1-badFileBeyondJustNE/float(goodFiles))

plt.show()



