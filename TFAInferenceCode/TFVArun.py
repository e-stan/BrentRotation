from TFVA import TFVA
import pickle
import numpy as np

#Get CS Matrix
file = open("TFA_CSMatrices/learnedCSval.csv","r")
delimiter = ","
CS = [[float(y) for y in x.split(delimiter)] for x in file.readlines()]

#Get TFA matrix
file = open("TFAFoundMatrix.csv","r")
TFAOrginial = [[float(y) for y in x.split(delimiter)] for x in file.readlines()]

#get GE Matrix
file = open("GERelevantHolstegeTrainedOnZev.csv","r")
GE = [[float(y) for y in x.split(delimiter)] for x in file.readlines()]

GES = np.matmul(CS,TFAOrginial)
print(np.mean(np.abs(np.subtract(GE,GES))))

print(TFVA(CS,GE,TFAOrginial))



