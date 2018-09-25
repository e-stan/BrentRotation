import numpy as np

file = open("learnedTFAvalQuantile5iteration100.csv","r")
delimiter = ","
matrix = [[float(y) for y in x.split(delimiter)] for x in file.readlines()]
print(len(matrix[0]))

columnLabelFile = "deleteome_all_mutants_controls.tsv"

columnLabels = open(columnLabelFile,"r").readline().split("\t")[3:]
print(len(columnLabels)/3.0)

realColumnLabels = []
perSample = 0
for column in columnLabels:
	if perSample == 2:
		try:
			column.index("vs. wt")
			realColumnLabels.append(str.lower(column.replace(" vs. wt","")[:-4])) 
		except:
			realColumnLabels.append(str.lower(column.replace(" vs wt","")[:-4]))
		perSample = 0 
	else:
		perSample += 1
print(len(realColumnLabels))

kinasesOfInterestFile = "finalKinasesToDetermineInteractions2.txt"

common2Systematic = {x.rstrip().split(" ")[0]:x.rstrip().split(" ")[-1] for x in open(kinasesOfInterestFile,"r").readlines()[1:]}
systematic2common = {value:key for key,value in common2Systematic.items()}

columnsOfInterest = []
columnsOfInterestNames = []
for x in range(len(realColumnLabels)):
	if realColumnLabels[x] in common2Systematic or realColumnLabels[x] in systematic2common:
		columnsOfInterest.append(x)
		if realColumnLabels[x] in common2Systematic: 
			columnsOfInterestNames.append(common2Systematic[realColumnLabels[x]])
		else:
			columnsOfInterestNames.append(realColumnLabels[x])

print(len(common2Systematic))
print(len(columnsOfInterest))
print(len(columnsOfInterestNames))
columnsOfInterestNames.append("wt")

wt = [np.mean([matrix[row][col] for col in [-1,-2,-3]]) for row in range(len(matrix))] #average the last 3 as there is a not a pure WT sample
matrix = [[matrix[row][col] for col in columnsOfInterest]+[wt[row]] for row in range(len(matrix))]

print(len(matrix[0]))


rowLabelFile = "TFListZev.csv"

rowLabels = [str.lower(x.rstrip()) for x in open(rowLabelFile,"r").readlines()]

print(len(rowLabels))


#Divide TF activities

matrix = [[matrix[row][col]/matrix[row][-1] for col in range(len(matrix[row]))] for row in range(len(matrix))]
print(len(matrix))
outputFile = open("NormalizedTFAMatrixOfInterest.csv","w")
outputFile.write("TF")
[outputFile.write(","+kinase) for kinase in columnsOfInterestNames]
for row in range(len(matrix)-1):
	outputFile.write("\n"+rowLabels[row])
	[outputFile.write(","+str(colEntry)) for colEntry in matrix[row]]



