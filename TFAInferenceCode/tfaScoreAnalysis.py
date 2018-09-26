import matplotlib.pyplot as plt

expectedResults = [x.rstrip().split()[:3] for x in open("../SelectedData/kinaseTFInteractionsWQuantititation.txt","r").readlines()[1:]]
print(expectedResults)
expectedResults = [[kinase+"_"+tf,score] for kinase,tf,score in expectedResults]
interactions = {x[0]:-1 for x in expectedResults}
for inter,score in expectedResults:
    if score == "BG:": score = 2.0
    interactions[inter] = score


foundResults = [x.rstrip().split(",") for x in open("NormalizedTFAMatrixOfInterest.csv","r").readlines()]
kinases = foundResults[0][1:]
print(kinases)
foundResults=foundResults[1:]
TFs =  [x[0] for x in foundResults]
print TFs
matrix = [[float(y) for y in x[1:]] for x in foundResults]
interactionsFound = {x[0]:-1 for x in expectedResults}
allInteractions = []
for tf in range(len(matrix)):
    for kinase in range(len(matrix[tf])):
        interactionsFound[kinases[kinase]+"_"+TFs[tf]] = matrix[tf][kinase]
        allInteractions.append(matrix[tf][kinase])

allInteractions.sort(reverse=True)
scoreCutoff = allInteractions[int(.05*len(allInteractions))]
print(scoreCutoff,int(.05*len(allInteractions)),len(allInteractions))
uniqueInteractions = [inter for inter in interactions]
interactionsCodeMap = {inter:code for inter,code in zip(uniqueInteractions,range(len(interactions)))}
print(len(uniqueInteractions),len([interactions[inter] for inter in uniqueInteractions if interactionsFound[inter] >= scoreCutoff]))
plt.scatter([interactions[inter] for inter in uniqueInteractions if interactionsFound[inter] >= scoreCutoff],[interactionsFound[inter] for inter in uniqueInteractions if interactionsFound[inter] >= scoreCutoff],c = "b")
plt.figure()
plt.hist(allInteractions, bins = 40)
plt.show()


