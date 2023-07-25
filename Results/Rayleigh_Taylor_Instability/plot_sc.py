import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import matplotlib.animation as animation
import sys
import glob


file_name=sys.argv[1]
file_name_split=file_name.split("/")
files=glob.glob(file_name_split[0]+"/"+file_name_split[1]+"/sc_file_*")
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
print('Results for Theta Evolution')
print('Reading mesh....')
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
print('Done Mesh Generation')
print('Reading Theta Results.....', end=" ");

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

print ('Done reading theta data.')	
times=np.array(times)

ind = np.argsort(times)
times=np.asarray([times[i] for i in ind]).flatten()
Data=[sc_Data[i] for i in ind]
print('Generating Animation...');

frames=[]
fig=plt.figure();



def animate(i):
	fig.clear()
	plt.suptitle("time=%3f"% (times[i]));
	ax_actual=fig.add_subplot(132)
	c=ax_actual.tripcolor(x_points, y_points,triangles, facecolors=np.asarray(Data[i]).T, edgecolors='none')
	cbar0=fig.colorbar(c,ax=ax_actual, orientation='vertical')
	#ax_actual.set_xlim([0,np.max(x_points)])
	ax_actual.set_ylim([0.5*np.max(y_points),np.max(y_points)])
	plt.gca().set_aspect('equal')
	plt.draw();

	
ani = animation.FuncAnimation(fig, animate,save_count=np.size(times))
#ani = animation.FuncAnimation(fig, animate,save_count=10)
writergif = animation.PillowWriter(fps=10) 
anim_name=file_name_split[0]+"/"+file_name_split[1]+"_theta_animation.gif"
ani.save(anim_name, writer='pillow')
print('Animation saved as: '+anim_name)
print('Done')
	
	


