import matplotlib.pyplot as plt
import copy


def flatten(l):
	return [item for sublist in l for item in sublist]

zevTFFile = "TFListZev.txt"
yeastKIDFile = "yeastKIDKinaseInteractions.txt"

bioGridFoundInteractions = [x.rstrip().split("\t")[:2]+x.rstrip().split("\t")[-3:-2] for x in open("finalBioGRIDYeastKinaseTFInteractions.tsv","r").readlines()]

kinaseTFInteractions = {}

for tf,kinase,score in bioGridFoundInteractions[1:]:
    if kinase not in kinaseTFInteractions: kinaseTFInteractions[kinase] = {}
    kinaseTFInteractions[kinase][tf] = "BG: "+score

kinaseTFInteractionsBIOGRIDOnly = copy.deepcopy(kinaseTFInteractions)

TFs = [str.lower(x.rstrip()) for x in open(zevTFFile,"r").readlines()][:-1]

file = open(yeastKIDFile,"r")
[file.readline() for x in range(8)]
headers = file.readline()
relevantInteractions = []
scores = []
for line in file:
    interaction = line.split("\t")
    if str.lower(interaction[2]) in TFs or str.lower(interaction[0]) in TFs:
        if str.lower(interaction[3]) not in kinaseTFInteractions: kinaseTFInteractions[str.lower(interaction[3])] = {}
        kinaseTFInteractions[str.lower(interaction[3])][str.lower(interaction[2])] = float(interaction[-1])
        relevantInteractions.append(line)

print("Number of interactions Total: ",sum([len(kinaseTFInteractions[kinase]) for kinase in kinaseTFInteractions]))

outputKinaseTFFile = open("kinaseTFInteractionsWQuantititation.txt","w")

unique2BioGrid = 0
unique2yeastKID = 0
scoresYeastKID = []
scoresShared = []
for kinase in kinaseTFInteractions:
    for tf in kinaseTFInteractions[kinase]:
        if str(kinaseTFInteractions[kinase][tf]).startswith("BG:"):
            unique2BioGrid+=1
            print(kinase,tf)
        elif kinase not in kinaseTFInteractionsBIOGRIDOnly or tf not in kinaseTFInteractionsBIOGRIDOnly[kinase]:
            unique2yeastKID+=1
            scoresYeastKID.append(kinaseTFInteractions[kinase][tf])
        else:
            scoresShared.append(kinaseTFInteractions[kinase][tf])
	outputKinaseTFFile.write(kinase+" "+tf+" "+str(kinaseTFInteractions[kinase][tf])+"\n")

print("Number of interactions unqiue to BioGRID: ",unique2BioGrid)
print("Number of interactions unique to yeastKID: ",unique2yeastKID)
print("Number of YeastKID interactions: ",len(relevantInteractions))
print("Number of BioGRID interactions: ",sum([len(kinaseTFInteractionsBIOGRIDOnly[kinase]) for kinase in kinaseTFInteractionsBIOGRIDOnly]))

relevantInteractions.sort(key=lambda x:x.split("\t")[-1])
outfile = open("relevantYeastKIDInteractions.txt","w")
outfile.write(headers)
[outfile.write(line) for line in relevantInteractions]


plt.hist(scoresYeastKID, color = "b")
plt.xlabel("Log Likelihood Score")
plt.ylabel("Frequency")
plt.title("YeastKID Only Kinase/Phosphatase TF Interactions Scores")
plt.figure()
plt.hist(scoresShared,color = "r")
plt.xlabel("Log Likelihood Score")
plt.ylabel("Frequency")
plt.title("BioGRID,YeastKID Shared Kinase/Phosphatase TF Interactions Scores")

plt.figure()
plt.hist([len(kinaseTFInteractions[kinase]) for kinase in kinaseTFInteractions])
plt.xlabel("Number of TFs Interacting with Kinase/Phosphatase")
plt.ylabel("Frequency")
plt.title("Distribution of Kinase/Phosphatase Interactions")

print("# of Kinases: ",len(kinaseTFInteractions))
print("# of TFs: ",len(set(flatten([[x for x in kinaseTFInteractions[kinase]] for kinase in kinaseTFInteractions]))))

plt.show()


