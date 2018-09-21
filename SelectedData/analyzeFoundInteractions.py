import matplotlib.pyplot as plt


def flatten(l):
	return [item for sublist in l for item in sublist]



kinaseListFile = "finalKinasesToDetermineInteractions.txt"
zevTFFile = "TFListZev.txt"

TFs = [str.lower(x.rstrip()) for x in open(zevTFFile,"r").readlines()][:-1]

kinases = [str.lower(x.split()[0]) for x in open(kinaseListFile,"r").readlines()[1:]]

filename = "bioGridRelevantInteractions.txt"
interactions = [x.rstrip().split("\t") for x in open(filename,"r")]
for x in range(len(interactions)): interactions[x][9:11] = [inter.split("|") for inter in interactions[x][9:11]]
headers = interactions[0]
interactions = interactions[1:]
kinaseTFInteractions = {tf:[] for tf in TFs}
interactionCount = 0
interactionRecord = []
for interaction in interactions:
	tf = [inter for inter in interaction[5:7] if str.lower(inter) in TFs]
	if len(tf) != 1: print("error")
	else: tf = str.lower(tf[0])
	kinase = [inter for inter in (interaction[5:9]+flatten(interaction[9:11])) if str.lower(inter) in kinases and inter != tf]
	if len(kinase) != 1: print("error",kinase,tf)
	elif interaction[12] == "physical":
		kinase = [str.lower(inter) for inter in interaction[5:7] if str.lower(inter) != tf][0]
		kinaseTFInteractions[tf].append(kinase)
		interactionCount+=1
		interactionRecord.append([tf,kinase,interaction[13],interaction[14],interaction[17],interaction[18],interaction[19],interaction[20]])
interactionRecord.sort(key=lambda x:x[0])
outputFileOfInteraction = open("finalBioGRIDYeastKinaseTFInteractions.tsv","w")
outputFileOfInteraction.write("TF\tKinase\tAuthor\tPubmed ID\tThroughput\tScore\tModification\tPhenotypes")
for entry in interactionRecord:
	outputFileOfInteraction.write("\n")
	for info in entry:
		outputFileOfInteraction.write(info+"\t")
		


print("Number of TFs with at least 1 relevant kinase/phosphatase interaction: ",len([inter for inter in kinaseTFInteractions if len(kinaseTFInteractions[inter]) > 0]))
print("Number of TFs in ZEV trained network: ", len(TFs))
interKinases = set()
for tf in kinaseTFInteractions:
	interKinases = interKinases.union(set(kinaseTFInteractions[tf]))
print("Number of Unique Kinases/Phosphatases Interacting: ",len(interKinases))
print("Number of Unique Kinases/Phosphatases in Kemmeren Dataset: ",len(kinases))

print("Total Number of Interactions: ",interactionCount)
plt.hist([len(kinaseTFInteractions[tf]) for tf in kinaseTFInteractions],bins = 35)
plt.xlabel("Number of Kinases/Phosphatases Interacting with TF")
plt.ylabel("Frequency")
plt.title("Distribution of TF Interactions")

kinaseInteractions = {}
for tf in kinaseTFInteractions:
	for kinase in kinaseTFInteractions[tf]:
		if kinase not in kinaseInteractions: kinaseInteractions[kinase] = []
		kinaseInteractions[kinase].append(tf)
plt.figure()
plt.hist([len(kinaseInteractions[kinase]) for kinase in kinaseInteractions],bins = 50)
plt.xlabel("Number of TFs Interacting with Kinase/Phosphatase")
plt.ylabel("Frequency")
plt.title("Distribution of Kinase/Phosphatase Interactions")



plt.show()
