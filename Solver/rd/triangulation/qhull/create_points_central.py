import numpy as np
from random import *

NX = 64
NY = 64

NTOT = 0

XMIN = -0.5
XMAX = 0.5

YMIN = -0.5
YMAX = 0.5

DX = (XMAX - XMIN)/float(NX)
DY = (YMAX - YMIN)/float(NY)

REF = 1

f = open('points.txt','w')

# f.write(str(2) + ' not_rbox ' + str(NTOT) + ' D2' +'\n')
# f.write(str(NTOT) + '\n')

for i in range(NX):
	X = float(i)/float(NX)*(XMAX - XMIN) + XMIN
	for j in range(NY):
		Y = float(j)/float(NY)*(YMAX - YMIN) + YMIN

		RANDX = 2.0*random()-1.0
		RANDY = 2.0*random()-1.0

		X = X + 0.05*RANDX
		Y = Y + 0.05*RANDY

		if X < -0.5:
			X = -0.5

		if X > 0.5:
			X = 0.5

		if Y < -0.5:
			Y = -0.5

		if Y > 0.5:
			Y = 0.5

		R = np.sqrt(X*X + Y*Y)

		if R < 0.2:
			for k in range(REF):
				RANDX = 2.0*random()-1.0
				RANDY = 2.0*random()-1.0

				XZOOM = X + DX*RANDX
				YZOOM = Y + DY*RANDY

				f.write(str(XZOOM)+"\t"+str(YZOOM)+'\n')
				NTOT = NTOT + 1


		f.write(str(X)+"\t"+str(Y)+'\n')
		NTOT = NTOT + 1

f.close()

print("2 create_points_central D2", NTOT)
print(NTOT)