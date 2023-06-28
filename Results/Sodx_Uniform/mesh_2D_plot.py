import numpy as np
import matplotlib.pyplot as plt
import sys

print('Creating Mesh Plot....')
case=sys.argv[1];
data = open(case+'/B_Scheme/active_0.txt')
lines=data.readlines()

data_mesh= open(case+'/Delaunay2D.txt')
xmin,ymin,xmax,ymax=data_mesh.readline().split()

#print(xmin,ymin,xmax,ymax)

triangles=[]
x_points=[]
x_points=np.array(x_points)
y_points=[]
y_points=np.array(y_points)
for i in range(np.size(lines)-1):
	X0,Y0,X1,Y1,X2,Y2,tbin,u00,u01,u02 = lines[i+1].split()
	x_array=np.array([float(X0),float(X1),float(X2)])
	y_array=np.array([float(Y0),float(Y1),float(Y2)])
	
	triangle=[];
	for j in range(3):
		values=np.where(x_points==x_array[j])[0]
		flag=True;
		for k in values:
			if (y_points[k]==y_array[j]): 
				triangle.append(k);
				flag=False
				break;
		if flag:
			x_points=np.append(x_points,x_array[j])	
			y_points=np.append(y_points,y_array[j])
			triangle.append(np.size(x_points)-1);	

	triangles.append(triangle)






fig, ax = plt.subplots()
ax.triplot(x_points, y_points,triangles,'.-',markersize=0.25,linewidth=0.5,label=((r'$n=%.0f$'%np.size(x_points))))
ax.set_xlim(float(xmin),float(xmax))
ax.set_ylim(float(ymin),float(ymax))
ax.set_xlabel(r'$x$')
ax.set_ylabel(r'$y$')
ax.legend(loc='upper right')
plt.savefig(case+'/'+case+'.svg')
print('Plot Saved at ' + case+'/'+case+'.svg')
print('Done.')
#plt.show()
