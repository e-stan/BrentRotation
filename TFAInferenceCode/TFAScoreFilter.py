import numpy as np
import pickle 

#file = open("learnedTFAvalQuantile5iteration100.csv","r")
file = open("TFA_CSMatrices/learnedTFAval.csv","r")
#file = open("learnedTFAvalSignConst.csv","r")
#file = open("TFVAResult.csv","r")

delimiter = ","
matrix = [[float(y) for y in x.split(delimiter)] for x in file.readlines()]
print(len(matrix[0]))

columnLabelFile = "deleteome_all_mutants_controls.tsv"

geneExpression = open(columnLabelFile,"r")

columnLabels = geneExpression.readline().split("\t")[3:]


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

wt = [np.mean([matrix[row][col] for col in [-1]]) for row in range(len(matrix))] #average the last 3 as there is a not a pure WT sample
matrix = [[matrix[row][col] for col in columnsOfInterest]+[wt[row]] for row in range(len(matrix)-1)]

matrix.append([1 for _ in range(len(matrix[0]))])

pickle.dump(matrix,open("TFAFound.pkl","wb"))
file = open("TFAFoundMatrix.csv","w")
file.write(str(matrix[0][0]))
[file.write(","+str(val)) for val in matrix[0][1:]]
for row in matrix[1:]:
	file.write("\n"+str(row[0]))
	for val in row[1:]:
		file.write(","+str(val))
file.close()
file = open("kinaseNamesForColumns","w")
[file.write(kinase+"\n") for kinase in columnsOfInterestNames]
file.close()
pickle.dump(columnsOfInterestNames,open("kinaseNamesOfInterest.pkl","wb"))

print(len(matrix[0]))

matrix = matrix[:-1]


rowLabelFile = "TFListZev.csv"

rowLabels = [str.lower(x.rstrip())[1:-1] for x in open(rowLabelFile,"r").readlines()]

print(len(rowLabels))


#Divide TF activities

def standardScore(matrix):
	return [[(matrix[row][col]-np.mean(matrix[row]))/np.std(matrix[row]) for col in range(len(matrix[row]))] for row in range(len(matrix))]

def absoluteFoldChange(matrix):
	return [[np.abs(matrix[row][col]/matrix[row][-1]-1) for col in range(len(matrix[row]))] for row in range(len(matrix))]

def rankingChange(matrix):
	tempMat = np.transpose(matrix)
	for x in range(len(tempMat)):
		temp = list(tempMat[x])
		temp = [abs(np.log2(temp[i]/matrix[i][-1])) for i in range(len(temp))]
		indices = range(len(temp))
		indices.sort(key=lambda x:temp[x],reverse=True)
		tempMat[x] = indices
	return np.transpose(tempMat).tolist()

def valueRelative2Max(matrix):
	tempMat = np.transpose(matrix)
	for x in range(len(tempMat)):
		temp = list(tempMat[x])
		temp = [(temp[i]/matrix[i][-1]) for i in range(len(temp))]
		tempMat[x] = [i/max(temp) for i in temp]
	return np.transpose(tempMat).tolist()

#FOR TFVA ONLY
file = open("TFVAResult.csv","r")

delimiter = ","
matrix = [[float(y) for y in x.split(delimiter)] for x in file.readlines()[:-1]]

matrix = valueRelative2Max(matrix)

print(len(matrix))
outputFile = open("NormalizedTFAMatrixOfInterest.csv","w")
outputFile.write("TF")
[outputFile.write(","+kinase) for kinase in columnsOfInterestNames]
for row in range(len(matrix)):
	outputFile.write('\n'+rowLabels[row])
	[outputFile.write(","+str(colEntry)) for colEntry in matrix[row]]





