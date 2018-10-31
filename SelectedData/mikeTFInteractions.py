
bioGridFile = "BIOGRID-ALL-3.4.164.tab2.txt"
bioGridFile = "BIOGRID-ORGANISM-Saccharomyces_cerevisiae_S288c-3.4.164.tab2.txt"
zevTFFile = "TFListZev.txt"
kinaseListFile = "YeastKinomeList.txt"
yeastKIDFile = "yeastKIDCompletePhysical.txt"

kinaseTypesGood = ["KinaseCatalytic","KinaseMetabolic/Lipid","PhosphataseCatalytic","PhosphataseMetabolic/Lipid"]

common2SystematicName = {str.lower(x.rstrip().split(",")[0]):str.lower(x.rstrip().split(",")[1]) for x in open("YeastCommonAndSystematicGeneNames.csv","r").readlines()[1:]}
systematic2Common = {value:key for key,value in common2SystematicName.items()}

foundUniqueInteractions = set()

def flatten(l):
	return [item for sublist in l for item in sublist]

TFs = ["mig1","mig2","eds1","rgt1"]

kinases = [str.lower(x.split()[0]) for x in open(kinaseListFile,"r").readlines()[1:]
	if "".join(x.rstrip().split()[1:]) in kinaseTypesGood]

tempKinases = []
tempTFs = []
for kinase in kinases:
	if kinase in common2SystematicName:
		tempKinases.append(common2SystematicName[kinase])
	elif kinase in systematic2Common:
		tempKinases.append(kinase)
for tf in TFs:
	if tf in common2SystematicName:
		tempTFs.append(common2SystematicName[tf])
	elif tf in systematic2Common:
		tempTFs.append(tf)

TFs = tempTFs
kinases = tempKinases


file = open(bioGridFile,"r")
outputInteractionFile = open("bioGridRelevantInteractionsMikeTFs.tsv","w")
i = 0
outputInteractionFile.write(file.readline())
relevantInteractions = []

for line in file:
	interactors = [str.lower(x) for x in line.split("\t")[5:7]]
	if any(tf in interactors for tf in TFs) and any(k in interactors for k in kinases) and line.split("\t")[12] == "physical":
		relevantInteractions.append(line)
		tf = [tf for tf in TFs if tf in interactors][0]
		kinase = [kinase for kinase in kinases if kinase in interactors][0]
		foundUniqueInteractions.add(str.upper(systematic2Common[str.lower(tf)])+"/"+str.upper(systematic2Common[str.lower(kinase)]))

relevantInteractions.sort(key = lambda x : x.split("\t")[5])

[outputInteractionFile.write(x) for x in relevantInteractions]

outputInteractionFile.close()


outputInteractionFile = open("yeastKIDRelevantInteractionsMikeTFs.tsv","w")

file = open(yeastKIDFile,"r")
[file.readline() for x in range(8)]
headers = file.readline()
outputInteractionFile.write(headers)
relevantInteractions = []
scores = []
for line in file:
    interaction = line.split("\t")
    if str.lower(interaction[2]) in TFs and str.lower(interaction[3]) in kinases:
        relevantInteractions.append(line)
        foundUniqueInteractions.add(str.upper(systematic2Common[str.lower(interaction[2])])+"/"+str.upper(systematic2Common[str.lower(interaction[3])]))

relevantInteractions.sort(key = lambda x : x.split("\t")[2])
[outputInteractionFile.write(x) for x in relevantInteractions]


foundUniqueInteractions = list(foundUniqueInteractions)
foundUniqueInteractions.sort()
for x in foundUniqueInteractions: print x

outputInteractionFile.close()