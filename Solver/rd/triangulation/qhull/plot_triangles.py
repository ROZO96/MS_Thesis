import numpy as np
import matplotlib.pyplot as plt

data_points = np.loadtxt('points.txt',skiprows=2)
data_triangles = np.loadtxt('rand_triangles.txt',skiprows=1)

x, y = data_points[:,0]+0.5, data_points[:,1]+0.5

n_vert, v0, v1, v2 = data_triangles[:,0], data_triangles[:,1], data_triangles[:,2], data_triangles[:,3]

# print(x)
# print(y)

# print(v0)
# print(v1)

plt.scatter(x,y)

for i in range(len(n_vert)):
    plt.plot((x[int(v0[i])],x[int(v1[i])]),(y[int(v0[i])],y[int(v1[i])]),color='red')
    plt.plot((x[int(v0[i])],x[int(v2[i])]),(y[int(v0[i])],y[int(v2[i])]),color='red')
    plt.plot((x[int(v1[i])],x[int(v2[i])]),(y[int(v1[i])],y[int(v2[i])]),color='red')
    if i % 100 == 0:
    	print(i,"/",len(n_vert))

# i=1

# plt.plot((x[int(v0[i])],x[int(v1[i])]),(y[int(v0[i])],y[int(v1[i])]),color='b')
# plt.plot((x[int(v0[i])],x[int(v2[i])]),(y[int(v0[i])],y[int(v2[i])]),color='b')
# plt.plot((x[int(v1[i])],x[int(v2[i])]),(y[int(v1[i])],y[int(v2[i])]),color='b')

plt.xlim(0.1,0.9)
plt.ylim(0.1,0.9)

plt.show()
