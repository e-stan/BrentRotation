import GEOparse
import pandas

filename = "deleteome_responsive_mutants_controls.txt"

df = pandas.read_csv(filename,sep="\t",low_memory=False)

kinaseProteomeListFile = "YeastKinomeList.txt"

kinaseTypeMap = {str.lower(kinase):kinaseType for kinase,kinaseType in [x.split()[:2] for x in open(kinaseProteomeListFile,"r").readlines()[1:]]}

columnsOfInterest = [column for column in list(df.columns.values) if any(kinase in str(column) for kinase in kinaseTypeMap)]

print(len(set(columnsOfInterest))/float(len(kinaseTypeMap)))  #94% of kinases/proteases in the kinome had impact of gene expression and was removed in kememmeron dataset








