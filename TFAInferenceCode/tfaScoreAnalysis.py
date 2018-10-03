import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats
import pickle

expectedResults = [x.rstrip().split()[:3] for x in open("../SelectedData/kinaseTFInteractionsWQuantititation.txt","r").readlines()[1:]]


geneExpression = pickle.load(open("geneExpressionForKinaseOfInterest.pkl","rb"))
geneExpression = {common2SystematicName[key]:value for key,value in geneExpression.items() if key in common2SystematicName}
for kinase in geneExpression:
	for type in geneExpression[kinase]:
		geneExpression[kinase][type] = {common2SystematicName[key]:value for key,value in geneExpression[kinase][type].items() if key in common2SystematicName}


expectedResults = [[kinase+"_"+tf,score] for kinase,tf,score in expectedResults]
interactions = {x[0]:-1 for x in expectedResults}
for inter,score in expectedResults:
    if score == "BG:": score = 2.0
    interactions[inter] = float(score)


foundResults = [x.rstrip().split(",") for x in open("NormalizedTFAMatrixOfInterest.csv","r").readlines()]
kinases = foundResults[0][1:]

foundResults=foundResults[1:]
TFs =  [x[0] for x in foundResults]

matrix = [[float(y) for y in x[1:]] for x in foundResults]
interactionsFound = {x[0]:-1 for x in expectedResults}
allInteractions = []
for tf in range(len(matrix)):
    for kinase in range(len(matrix[tf])):
        interactionsFound[kinases[kinase]+"_"+TFs[tf]] = matrix[tf][kinase]
        allInteractions.append(matrix[tf][kinase])
print([inter for inter in interactionsFound if interactionsFound[inter] == -1])
allInteractions.sort(reverse=True)
alpha = .05 # for standard score
#alpha = .1 # for absolute fold change
scoreCutoffHigh = allInteractions[int(alpha/2*len(allInteractions))]
allInteractions.sort()
scoreCutoffLow = allInteractions[int(alpha/2*len(allInteractions))]
#scoreCutoffLow = -.01 # for absolute fold change
alpha=.05
print(scoreCutoffHigh,scoreCutoffLow,alpha,int(alpha*len(allInteractions)),len(allInteractions))

uniqueInteractions = [inter for inter in interactions]
#interactionsCodeMap = {inter:code for inter,code in zip(uniqueInteractions,range(len(interactions)))}


significantFoundInteractions = [inter for inter in uniqueInteractions if interactionsFound[inter] >= scoreCutoffHigh or interactionsFound[inter] <= scoreCutoffLow]
insignificantInteractions= [inter for inter in uniqueInteractions if inter not in significantFoundInteractions]



print(len(uniqueInteractions),len(significantFoundInteractions))
plt.scatter([interactions[inter] for inter in significantFoundInteractions],[interactionsFound[inter] for inter in significantFoundInteractions],c = "b",label = "Significant")
plt.scatter([interactions[inter] for inter in insignificantInteractions],[interactionsFound[inter] for inter in insignificantInteractions],c="r",label = "Insignficiant")
plt.legend()
expectedInteractionLogScores = [interactions[inter] for inter in uniqueInteractions]
minLogScore = int(min(expectedInteractionLogScores))
maxLogScore = int(max(expectedInteractionLogScores)+1)

plt.plot([x for x in range(minLogScore,maxLogScore)],[scoreCutoffHigh for _ in range(minLogScore,maxLogScore)],c="black")
plt.plot([x for x in range(minLogScore,maxLogScore)],[scoreCutoffLow for _ in range(minLogScore,maxLogScore)],c="black")

plt.xlabel("Experimentally Validated Log Likelihood Score")
plt.ylabel("Standard Activity Score ((x-Mean[X])/Std[X])") #for standard score
#plt.ylabel("Absolute Fold Change Relative to wt") #for absolute fold change


plt.figure()
plt.hist(allInteractions, bins = 40)
plt.xlabel("Score")
plt.ylabel("Frequency")

a = len(significantFoundInteractions)
b = len(uniqueInteractions)-a
c = int(alpha*len(allInteractions)) - a
d = len(allInteractions) - a - b - c

print(stats.fisher_exact([[a,b],[c,d]],alternative = "greater"))

GEImpact = {inter:-1 for inter in uniqueInteractions}
xvals = []
yvals = []
for inter in uniqueInteractions:
	kinase = inter.split("_")[0]
	allScoresForKinase = [geneExpression[kinase2]["expression"][kinase] for kinase2 in geneExpression]
	kianseStdScore = (geneExpression[kinase]["expression"][kinase] - np.mean(allScoresForKinase))/np.std(allScoresForKinase)
#	kinaseStdScore = abs(2**geneExpression[kinase]["ratio"][kinase]-1) #for absolute fold change only
	GEImpact[inter] = kianseStdScore
	xvals.append(abs(kianseStdScore))
	yvals.append(abs(interactionsFound[inter]))

plt.figure()
plt.scatter([abs(GEImpact[inter]) for inter in significantFoundInteractions],[abs(interactionsFound[inter]) for inter in significantFoundInteractions],c = "b",label = "Significant")
plt.scatter([abs(GEImpact[inter]) for inter in insignificantInteractions],[abs(interactionsFound[inter]) for inter in insignificantInteractions],c="r",label = "Insignficiant")
plt.legend()
plt.xlabel("Absolute Value of Standard Score for Kinase/Phosphotase Expression")
plt.ylabel("Absolute Value of Standard Score for TFA")
#plt.xlabel("Absolute Fold Change of Kinase/Phosphotase Expression vs WT") #only for absolute FC
#plt.ylabel("Absolute Fold Change of TFA vs WT")

positions = range(len(xvals))
positions.sort(key=lambda x: xvals[x])
print(positions)
xvals = [xvals[x] for x in positions]
yvals = [yvals[x] for x in positions]
slope,intercept,r,pvalue,_ = stats.linregress(xvals, yvals)
plt.plot(xvals,[x*slope+intercept for x in xvals],"r-")
print(r**2,pvalue)

plt.show()


