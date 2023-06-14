import numpy as np

points = np.loadtxt('points.txt',skiprows=2)
triangles = np.loadtxt('triangles.txt',skiprows=1)

x, y = points[:,0], points[:,1]
n_vert, vert0, vert1, vert2 = triangles[:,0], triangles[:,1], triangles[:,2], triangles[:,3]

t = open('triangles.txt','r')

n_triangles = t.readline()

f = open('ordered_triangles.txt','w')

f.write(n_triangles)

for i in range(len(n_vert)):
	x0 = x[int(vert0[i])]
	y0 = y[int(vert0[i])]
	x1 = x[int(vert1[i])]
	y1 = y[int(vert1[i])]
	x2 = x[int(vert2[i])]
	y2 = y[int(vert2[i])]

	l1x = x1 - x0
	l1y = y1 - y0

	l2x = x2 - x0
	l2y = y2 - y0

	cross = l1x*l2y - l1y*l2x

	if cross < 0.0:
		vert1[i] = vert1[i] + vert2[i]
		vert2[i] = vert1[i] - vert2[i]
		vert1[i] = vert1[i] - vert2[i]
		
	f.write(str(int(n_vert[i])) + "\t" + str(int(vert0[i])) + "\t" + str(int(vert1[i])) + "\t" + str(int(vert2[i])) + "\n")

f.close()