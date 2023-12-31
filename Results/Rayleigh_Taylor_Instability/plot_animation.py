import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import matplotlib.animation as animation
import sys
import glob
import os

file_name=sys.argv[1]
file_name_split=file_name.split("/")
files=sorted(glob.glob(file_name_split[0]+"/"+file_name_split[1]+"/snapshot_*"),key=os.path.getmtime)

fig = plt.figure()


def animate(i):
	fp=open(files[i])
	lines=fp.readlines()

	fig.clear()
	n_point,time=int(lines[0].split()[0]),float(lines[0].split()[1])
	#plt.title('Time= %f'%(time))
	plt.axis('off')
	data = []
	for string in lines[1:]:
    		data.append(string.split())
	data=np.asarray(data).astype(float)
	x,y=data[:,0],data[:,1]
	
	titles=['Density','Velocity x','Velocity v_y','Pressure P']
	
	for i in range(3):
		
		z=data[:,i]
	
		num_plot=219+i;
		#ax_actual=fig.add_subplot(num_plot)
		#ax_actual=fig.add_subplot(111)
		#ax_actual.set_title(titles[i-2])
		
		if (i==2):
			ax_actual=fig.add_subplot(111)
			ax_actual.set_aspect('equal')
			ax_actual.axis('equal')
			ax_actual.set_title('Time= %f'%(time))
			c=ax_actual.tricontourf(x,y,z,np.linspace(1, 3, 32),
                   	extend='both')
			cbar0=fig.colorbar(c,ax=ax_actual, orientation='vertical')
			cbar0.set_ticks(np.arange(1, 3.25, .25))
			ax_actual.set_xlim(x.min(),x.max())
			#ax_actual.set_ylim(y.min(),0.5*y.max())
			ax_actual.set_ylim(0.5*y.max(),y.max())
			fig.tight_layout()
			plt.gca().set_aspect('equal')
			plt.draw()
		'''
		else:
			c=ax_actual.tricontourf(x,y,z,levels=32)
			cbar0=fig.colorbar(c,ax=ax_actual, orientation='vertical')
		'''
		

	
	

ani = animation.FuncAnimation(fig, animate,save_count=np.size(files))
writergif = animation.PillowWriter(fps=10) 
anim_name=file_name_split[0]+"/"+file_name_split[1]+"_plot_animation_45.gif"
ani.save(anim_name, writer='pillow')
print('Animation saved as: '+anim_name)
print('Done')

