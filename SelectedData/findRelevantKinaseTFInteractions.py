
bioGridFile = "BIOGRID-ALL-3.4.164.tab2.txt"
bioGridFile = "BIOGRID-ORGANISM-Saccharomyces_cerevisiae_S288c-3.4.164.tab2.txt"
zevTFFile = "TFListZev.txt"
kinaseListFile = "finalKinasesToDetermineInteractions.txt"


def flatten(l):
	return [item for sublist in l for item in sublist]

TFs = [str.lower(x.rstrip()) for x in open(zevTFFile,"r").readlines()][:-1]

kinases = [str.lower(x.split()[0]) for x in open(kinaseListFile,"r").readlines()[1:]]

file = open(bioGridFile,"r")
outputInteractionFile = open("bioGridRelevantInteractions.txt","w")
i = 0
outputInteractionFile.write(file.readline())
relevantInteractions = []
for line in file:
	interactors = flatten([[str.lower(gene) for gene in x.split("|")] for x in (line.split("\t")[7:11]+line.split("\t")[5:7])])
	if any(tf in interactors for tf in TFs) and any(k in interactors for k in kinases):

		relevantInteractions.append(line)


relevantInteractions.sort(key = lambda x : x.split("\t")[5])

[outputInteractionFile.write(x) for x in relevantInteractions]


	
	



