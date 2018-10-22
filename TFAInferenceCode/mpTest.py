from TFHelp import *
import numpy as np

GE = [[1,2,3],[1,2,0],[2,2,2]]
CS = [[1,0,0],[0,1,0],[0,0,1]]
TFA = np.dot(np.linalg.pinv(CS),GE).tolist()



print(TFAHelper.TFVAMPRun(GE,CS,TFA,write=False,numCores = 3))


print(TFAHelper.TFVA(GE,CS,TFA,write=False))