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
print(triangles)


fig, ax = plt.subplots()
#triang=mtri.Triangulation(val_x, val_y, triangles)
triang= mtri.Triangulation(val_x, val_y)
print(triang.triangles)
ax.triplot(val_x, val_y,triangles,'o-',markersize=0.1,linewidth=0.5)
ax.triplot(triang,'o-',markersize=1,linewidth=0.5)
"""
for i in range(n_triangles):
    ax.plot((val_x[int(triangles[i,0])],val_x[int(triangles[i,1])]),(val_y[int(triangles[i,0])],val_y[int(triangles[i,1])]),color='red')
    ax.plot((val_x[int(triangles[i,0])],val_x[int(triangles[i,2])]),(val_y[int(triangles[i,0])],val_y[int(triangles[i,2])]),color='red')
    ax.plot((val_x[int(triangles[i,1])],val_x[int(triangles[i,2])]),(val_xy[int(triangles[i,1])],val_y[int(triangles[i,2])]),color='red')
"""
ax.set_xlim(float(xmin),float(xmax))
ax.set_ylim(float(ymin),float(ymax))

plt.show()
