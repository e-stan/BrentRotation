import pandas
import matplotlib.pyplot as plt
import numpy
import pickle

def flatten(l):
	return [item for sublist in l for item in sublist]


#read in data

filename = "deleteome_all_mutants_controls.tsv"

df = pandas.read_csv(filename,sep="\t",low_memory=False)

df["geneSymbol"] = df["geneSymbol"].str.lower()

#get list of kinase/phosphotases in yeast

kinaseProteomeListFile = "finalKinasesToDetermineInteractions2.txt"

kinases = [[x.rstrip().split()[0],x.rstrip().split()[-1]] for x in open(kinaseProteomeListFile,"r").readlines()[1:]]
common2SystematicKinase = {key:value for key,value in kinases}
kinases = flatten(kinases)

#determine which kinases/phosphatases were knocked out in experiment
columnsTotal = list(df.columns.values)[3:]
for column,line in zip(list(df.columns.values)[3:],range(len(list(df.columns.values)))):
	try:
		column.index("vs. wt")
	except:
		columnsTotal[line] = columnsTotal[line].replace("vs wt","vs. wt")

columnsOfInterestFullName = [column for column in columnsTotal if any(kinase in str(column) for kinase in kinases)]
columnsOfInterest = [column[0:column.index("-del vs. wt")] for column in columnsTotal if any(kinase in str(column) for kinase in kinases)]

columnsOfInterestFullName = [columnsOfInterestFullName[x] for x in range(len(columnsOfInterestFullName)) if any(columnsOfInterest[x] == kinase for kinase in kinases)]
columnsOfInterest = [columnsOfInterest[x] for x in range(len(columnsOfInterest)) if any(columnsOfInterest[x] == kinase for kinase in kinases)]

df['systematicName'] = df['systematicName'].str.lower()


df = df.set_index('systematicName')

df = df.rename({numpy.nan:"dtype"},axis="index")
df2 = df.drop(["dtype"])

common2SystematicName = {str.lower(x.rstrip().split(",")[0]):str.lower(x.rstrip().split(",")[1]) for x in open("YeastCommonAndSystematicGeneNames.csv","r").readlines()[1:]}
systematic2Common = {value:key for key,value in common2SystematicName.items()}

dataOfInterest = {}
for column,gene in zip(columnsOfInterestFullName,columnsOfInterest):
	if gene in common2SystematicKinase: gene = common2SystematicKinase[gene]
	if gene not in dataOfInterest: dataOfInterest[gene] = {}	
	if str(df.at["dtype",column]) == "p_value": dataOfInterest[gene]["p-value"] = {key:float(value) for key,value in df2[column].to_dict().items()}
	elif str(df.at["dtype",column]) == "M": dataOfInterest[gene]["ratio"] = {key:float(value) for key,value in df2[column].to_dict().items()} #log2(mutant/WT)
	elif str(df.at["dtype",column]) == "A": dataOfInterest[gene]["expression"] = {key:float(value) for key,value in df2[column].to_dict().items()} #log2(wt)

temp = common2SystematicKinase[columnsOfInterest[0]]
dataOfInterest["wt"] = {}
dataOfInterest["wt"]["expression"] = {gene:dataOfInterest[temp]["expression"][gene]/(2**dataOfInterest[temp]["ratio"][gene]) for gene in dataOfInterest[temp]["ratio"]}

allData = {}
for column in columnsTotal:
	if column not in allData: allData[column] = {}
	if str(df.at["dtype",column]) == "p_value": allData[gene]["p-value"] = {key:float(value) for key,value in df2[column].to_dict().items()}
	elif str(df.at["dtype",column]) == "M": allData[gene]["ratio"] = {key:float(value) for key,value in df2[column].to_dict().items()} #log2(mutant/WT)
	elif str(df.at["dtype",column]) == "A": allData[gene]["expression"] = {key:float(value) for key,value in df2[column].to_dict().items()} #log2(wt)

temp = columnsTotal[-1]

allData[temp]["expression"] = {gene:allData[temp]["expression"][gene]/(2**allData[temp]["ratio"][gene]) for gene in allData[temp]["ratio"]}

allData = {kinase:{key:value for key,value in allData[kinase]["expression"].items()} for kinase in allData}



dataOfInterest = {kinase:{key:value for key,value in dataOfInterest[kinase]["expression"].items()} for kinase in dataOfInterest}

kinaseNames = pickle.load(open("kinaseNamesOfInterest.pkl","rb"))
geneNames = [str.lower(x.strip())[1:-1] for x in open("geneListZev.csv","r").readlines()]

GEMatrix = [[] for _ in geneNames]

for x in range(len(geneNames)):
	GEMatrix[x] = [dataOfInterest[kinase][geneNames[x]] for kinase in kinaseNames]


file = open("GERelevantHolstegeTrainedOnZev.csv","w")

file.write(str(GEMatrix[0][0]))
[file.write(","+str(val)) for val in GEMatrix[0][1:]]
for row in GEMatrix[1:]:
	file.write("\n"+str(row[0]))
	for val in row[1:]:
		file.write(","+str(val))

GEMatrix = [[] for _ in geneNames]

for x in range(len(geneNames)):
	GEMatrix[x] = [allData[kinase][geneNames[x]] for kinase in columnsTotal]


file = open("GEAll.csv","w")

file.write(str(GEMatrix[0][0]))
[file.write(","+str(val)) for val in GEMatrix[0][1:]]
for row in GEMatrix[1:]:
	file.write("\n"+str(row[0]))
	for val in row[1:]:
		file.write(","+str(val))



