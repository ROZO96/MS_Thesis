import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.ticker as ticker
import sys
import glob
import os

file_name=sys.argv[1]
file_name_split=file_name.split("/")
files=sorted(glob.glob(file_name_split[0]+"/"+file_name_split[1]+"/snapshot_*"),key=os.path.getmtime)

fig=plt.figure()


def animate(i):
	print('Creating frame %d of %d ...'%(i+1,np.size(files)),end=' ') 
	fp=open(files[i])
	lines=fp.readlines()
	n_point,time=int(lines[0].split()[0]),float(lines[0].split()[1])

	data = []
	for string in lines[1:]:
		data.append(string.split())
	data=np.asarray(data).astype(float)

	x,y=data[:,0],data[:,1]
	if 'vel' in file_name_split[1].split("_"):
		y=-y+2
	rho=data[:,2]

	ind_2=np.where((y>0.5*y.max()));

	x=np.asarray([x[i] for i in ind_2]).flatten()
	y=np.asarray([y[i] for i in ind_2]).flatten()
	rho=np.asarray([rho[i] for i in ind_2]).flatten()
	
	fig.clear()
	fig.suptitle('time=%1.3f'%(time))
	ax_0=fig.add_subplot(121)
	c=ax_0.tricontour(x,y,rho,levels=[(rho.min()+rho.max())/2.0],colors='lightseagreen') 
	#ax_0.clabel(c, inline=True, fontsize=10,colors='lightseagreen')
	ax_0.set_xlim(x.min(),x.max())
	ax_0.xaxis.set_major_locator(ticker.NullLocator())
	ax_0.yaxis.set_major_locator(ticker.NullLocator())
	ax_0.set_aspect('equal')

	ax_1=fig.add_subplot(122)
	c=ax_1.tricontourf(x,y,rho,levels=32) 
	plt.colorbar(c)
	ax_1.set_xlim(x.min(),x.max())
	ax_1.xaxis.set_major_locator(ticker.NullLocator())
	ax_1.yaxis.set_major_locator(ticker.NullLocator())
	ax_1.set_aspect('equal')
	plt.draw()
	print('Done') 

ani = animation.FuncAnimation(fig, animate,save_count=np.size(files))
#ani = animation.FuncAnimation(fig, animate,save_count=10)
writergif = animation.PillowWriter(fps=10) 
anim_name=file_name_split[0]+"/"+file_name_split[1]+"_density_plot_animation.gif"
ani.save(anim_name, writer='pillow')
print('Animation saved as: '+anim_name)
print('Done')


