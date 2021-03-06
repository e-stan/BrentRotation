import pandas
import matplotlib.pyplot as plt
import numpy
import pickle

#read in data

filename = "deleteome_responsive_mutants_controls.txt"
#filename = "deleteome_all_mutants_controls.tsv"
common2Syst = {str.lower(x.rstrip().split(",")[0]):str.lower(x.rstrip().split(",")[1]) for x in open("../TFAInferenceCode/YeastCommonAndSystematicGeneNames.csv","r").readlines()[1:]}
syst2common = {value:key for key,value in common2Syst.items()}


df = pandas.read_csv(filename,sep="\t",low_memory=False)


df["geneSymbol"] = df["geneSymbol"].str.lower()


#get list of kinase/phosphotases in yeast

kinaseProteomeListFile = "YeastKinomeList.txt"
def mergeListOfStrings(l,sep = " "):
	string = ""
	for st in l: 
		string+=(sep+st)
	return string

kinaseTypeMap = {str.lower(info[0]):mergeListOfStrings(info[1:]) for info in [x.split() for x in open(kinaseProteomeListFile,"r").readlines()[1:]]}

#determine which kinases/phosphatases were knocked out in experiment
columnsTotal = list(df.columns.values)[3:]
for column,line in zip(list(df.columns.values)[3:],range(len(list(df.columns.values)))):
	try:
		column.index("vs. wt")
	except:
		columnsTotal[line] = columnsTotal[line].replace("vs wt","vs. wt")
columnsOfInterestFullName = [column for column in columnsTotal if any(kinase in str(column) for kinase in kinaseTypeMap)]
columnsOfInterest = [column[0:column.index("-del vs. wt")] for column in columnsTotal if any(kinase in str(column) for kinase in kinaseTypeMap)]
columnsTotal = [column[0:column.index("vs. wt")-1] for column in columnsTotal]
print(len(set(columnsOfInterest))/float(len(kinaseTypeMap)))  #31% of kinases/proteases in the kinome had impact of gene expression and was removed in kememmeron dataset
print(len(set(columnsOfInterest)))# (81 kinases/proteases)
print(len(set(list(columnsTotal)))) #this is slightly off from 700 due to the following  (703)


#Confirm knockdown of selected gene

df = df.set_index("geneSymbol")

df = df.rename({numpy.nan:"dtype"},axis="index")

df2 = df.drop(["dtype"])
dataOfInterest = {}
for column,gene in zip(columnsOfInterestFullName,columnsOfInterest):
	if gene in list(df2.index.values):
		if gene not in dataOfInterest: dataOfInterest[gene] = {}	
		if str(df.at["dtype",column]) == "p_value": dataOfInterest[gene]["p-value"] = {key:float(value) for key,value in df2[column].to_dict().items()}
		elif str(df.at["dtype",column]) == "M": dataOfInterest[gene]["ratio"] = {key:float(value) for key,value in df2[column].to_dict().items()} #log2(mutant/WT)
		elif str(df.at["dtype",column]) == "A": dataOfInterest[gene]["expression"] = {key:float(value) for key,value in df2[column].to_dict().items()} #log2(wt)

dataOfInterest = {key:value for key,value in dataOfInterest.items() if key in kinaseTypeMap}
pickle.dump(dataOfInterest,open("geneExpressionForKinaseOfInterest.pkl","wb"))
plt.hist([dataOfInterest[gene]["ratio"][gene] for gene in dataOfInterest],bins = 25)
plt.xlabel("M (Log2(Mutant/WT)) Ratio For KO Gene")
plt.ylabel("Frequency")
plt.title("Kemmeren Kinase/Phosphotase KO Success")

#Quantify distribution of altered genes
plt.figure()
plt.hist([len([1 for gene2 in dataOfInterest[gene]["ratio"] if abs(dataOfInterest[gene]["ratio"][gene2])>numpy.log2(1.7) and dataOfInterest[gene]["p-value"][gene2] < .05])/float(len(dataOfInterest[gene]["ratio"])) for gene in dataOfInterest],bins = 100)

factoredGenes = []
for gene in open("geneListZev.csv","r").readlines():
	try:
		factoredGenes.append(str.lower(syst2common[str.lower(gene.rstrip())[1:-1]]))
	except:
		print gene 
print factoredGenes
# numDifferentiallyExpressedGenes = [[len([1 for gene2 in dataOfInterest[gene]["ratio"] if
# 			abs(dataOfInterest[gene]["ratio"][gene2])>numpy.log2(1.7) and gene2 in factoredGenes]),gene] for gene in dataOfInterest]
numDifferentiallyExpressedGenes = [[sum([abs(dataOfInterest[gene]["ratio"][gene2]) for gene2 in dataOfInterest[gene]["ratio"] if
			gene2 in factoredGenes]),gene] for gene in dataOfInterest]
numDifferentiallyExpressedGenes.sort(key = lambda x : x[0],reverse = True)
outfile = open("kinasesDifferentiallyExpressedGenes.txt","w")
for num,gene in numDifferentiallyExpressedGenes:
	outfile.write(str(num)+" "+gene+"\n")
outfile.close()
plt.xlabel("% of Genes Differentially Expressed (|M| > log2(1.7), p<.05)")
plt.ylabel("Frequency")
plt.title("Kemmeren Kinase/Phosphotase KO Differential Gene Expression")

#Select only kinases/phosphate K0 for which the gene is differentially expressed
qualityKinases2analyze = [gene for gene in dataOfInterest if abs(dataOfInterest[gene]["ratio"][gene]) > numpy.log2(1.7) and dataOfInterest[gene]["p-value"][gene] < .05]
print(len(qualityKinases2analyze)) #75 kinases/phosphotases total that meet QC standards


file=open("qualityKinaseKO.txt","w")
[file.write(gene+" "+kinaseTypeMap[gene]+"\n") for gene in qualityKinases2analyze]
file.close()

#show distribution of different types

differentKinTypesPresent = list(set([kinaseTypeMap[gene] for gene in qualityKinases2analyze]))
wedgeSizes = [len([1 for gene in qualityKinases2analyze if kinaseTypeMap[gene] == tp]) for tp in differentKinTypesPresent]
plt.figure()
plt.pie(wedgeSizes,labels = differentKinTypesPresent)


#Hand curated results

#file = open("qualityKinaseKOHandAnnotation.txt","r")
file = open("qualityKinaseKO.txt","r")

kinaseQualityHand = {str.lower(info[0]):mergeListOfStrings(info[1:]) for info in [x.split() for x in file.readlines()]}
kinaseTypesGood = [" Kinase Catalytic"," Kinase Metabolic/Lipid"," Phosphatase Catalytic"," Phosphatase Metabolic/Lipid"]
kinaseQualityHand = {key:value for key,value in kinaseQualityHand.items() if value in kinaseTypesGood}
differentKinTypesPresent = kinaseTypesGood
wedgeSizes = [len([1 for gene in kinaseQualityHand if kinaseQualityHand[gene] == tp]) for tp in differentKinTypesPresent]
print(len(kinaseQualityHand))
plt.figure()
plt.pie(wedgeSizes,labels = differentKinTypesPresent)

#Get systematic names


file.close()

file = open("finalKinasesToDetermineInteractions2.txt","w")
file.write("gene type\n")
for gene in kinaseQualityHand: 
	try:
		file.write(gene+" "+kinaseQualityHand[gene]+" "+common2Syst[gene]+"\n")
	except:
		if gene in syst2common: 
			file.write(syst2common[gene]+" "+kinaseQualityHand[gene]+" "+gene+"\n")
		else:
			print gene

file.close()

#plt.show()

