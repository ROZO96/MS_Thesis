import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri

data = open('Delaunay2D.txt')
lines=data.readlines()
xmin,ymin,xmax,ymax=lines[0].split()

print(xmin,ymin,xmax,ymax)

n_points=int(lines[2])
print(n_points)

n_triangles=int(lines[n_points+4])
print(n_triangles)

points = []
for string in lines[3:n_points+3]:
    points.append(string.split())
points=np.asarray(points)
val_x,val_y=points[:,0].astype(float),points[:,1].astype(float)

triangles = []
for string in lines[n_points+5:n_points+5+n_triangles]:
    triangles.append(string.split())
triangles=np.asarray(triangles).astype(int)


fig, ax = plt.subplots()
triang=mtri.Triangulation(val_x, val_y, triangles)
plt.triplot(triang,'o-',markersize=1,linewidth=0.1)
plt.xlim(xmin,xmax)
plt.ylim(ymin,ymax)

plt.show()
