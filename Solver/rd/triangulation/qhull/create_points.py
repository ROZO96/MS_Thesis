import numpy as np

NX = 64
NY = 64

NTOT = NX*NY

XMIN = -0.5
XMAX = 0.5

YMIN = -0.5
YMAX = 0.5

f = open('points.txt','w')

f.write(str(2) + ' not_rbox ' + str(NTOT) + ' D2' +'\n')
f.write(str(NTOT) + '\n')

for i in range(NX):
	for j in range(NY):
		if j % 2 == 0:
			X = float(i)/float(NX)*(XMAX - XMIN) + XMIN
		else:
			X = (float(i))/float(NX)*(XMAX - XMIN) + XMIN
		Y = float(j)/float(NY)*(YMAX - YMIN) + YMIN
		f.write(str(X)+"\t"+str(Y)+'\n')

f.close()
