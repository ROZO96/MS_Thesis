import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import matplotlib.animation as animation
import sys
import glob



files=glob.glob("Case_1/Bx_Scheme/sc_file_*")
fp=open(files[0])
lines=fp.readlines()
sc_Data=[];
times=[];
n_triangle,time=int(lines[0].split()[0]),float(lines[0].split()[1])
triangles=[]
x_points=[]
x_points=np.array(x_points)
y_points=[]
y_points=np.array(y_points)
sc_values=[]
for i in range(np.size(lines)-1):
	X0,Y0,X1,Y1,X2,Y2,SC = lines[i+1].split()
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
	sc_values.append(float(SC))
	
times.append(float(time))
sc_Data.append(sc_values)


for file_name in files[1:]:
	fp=open(file_name)
	lines=fp.readlines()
	n_triangle,time=int(lines[0].split()[0]),float(lines[0].split()[1])
	sc_current=[]
	for i in range(np.size(lines)-1):
		X0,Y0,X1,Y1,X2,Y2,SC = lines[i+1].split()
		sc_current.append(float(SC))	
	sc_Data.append(sc_current)
	times.append(float(time))

	
times=np.array(times)

ind = np.argsort(times)
times=np.asarray([times[i] for i in ind]).flatten()
Data=[sc_Data[i] for i in ind]
frames=[]
fig,ax=plt.subplots();



def animate(i):
	plt.clf() 
	plt.tripcolor(x_points, y_points,triangles, facecolors=np.asarray(Data[i]).T, edgecolors='k')
	plt.colorbar()
	plt.title("time=%f"% (times[i]));
	plt.xlim([0,2])
	plt.ylim([0,2])
	plt.draw();

	
	
ani = animation.FuncAnimation(fig, animate,save_count=np.size(times))
writergif = animation.PillowWriter(fps=30) 
ani.save('myAnimation.gif', writer='pillow')
##frame_one.save("my_awesome.gif", format="GIF", append_images=frames,save_all=True, duration=100, loop=0)
	
	


