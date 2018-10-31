import cobra

def getModel(filename):
    model = cobra.io.read_sbml_model(filename)

    # Correct metabolite ids:
    for met in model.metabolites:
        met.id = met.id.replace('__91__', '_')
        met.id = met.id.replace('__93__', '')

    model.solver = 'gurobi'
    return model



fluxDataFile = "exampleTestFluxData.csv"
delimiter = ","
GEMFileName = ""
error = .1

fluxData = [x.rstrip().split(delimiter) for x in open(fluxDataFile,"r").readlines()[1:]]

constrainingReactions = ['r_1761','r_1714','r_1672','r_1808'] #ethanol,glucose,CO2,glycerol

#GEMWithMedia = "yeastMinimalMediaGEM.xml"
GEMWithMedia = "yeastGEM.xml"


for sample in fluxData:
	model = getModel(GEMWithMedia)
	for reaction,flux in zip(constrainingReactions,[float(x) for x in sample[2:]]):
		model.reactions.get_by_id(reaction).upper_bound = flux +error
		model.reactions.get_by_id(reaction).lower_bound = flux - error
	
	sol = model.optimize()

	if sol == "optimal": cobra.io.write_sbml_model(model,sample[0]+"_"+sample[1]+".xml")
	else: print("error: ",sample[0]+"_"+sample[1]) 

