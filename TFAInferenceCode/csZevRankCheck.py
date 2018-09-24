import numpy
def setTo1or0(val):
	if val == 0:return 0
	else:return 1

matrix = [[setTo1or0(float(y)) for y in x.rstrip().split(",")] for x in open("signedBinaryCSZev.csv","r").readlines()]


print len(matrix[0])

print numpy.linalg.matrix_rank(matrix)
