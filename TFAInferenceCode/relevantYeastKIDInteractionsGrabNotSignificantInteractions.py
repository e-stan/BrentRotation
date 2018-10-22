data = [x.split(",") for x in open("../SelectedData/relevantYeastKIDInteractionsEthanAnnotation.csv","r").readlines()[1:]]
outfile = open("notFoundInteractionsFormatted.csv","w")
for x in data:
	for i in range(4):
		outfile.write(x[i]+",")
	for i in xrange(4,len(x)):
		if x[i] != "":
			outfile.write(" "+x[i])
