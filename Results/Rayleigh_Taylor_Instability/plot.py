import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate
import sys
import matplotlib as mpl

file_name=sys.argv[1]
fp=open(file_name)
lines=fp.readlines()

n_point,time=int(lines[0].split()[0]),float(lines[0].split()[1])
print(n_point,time)



data = []
for string in lines[1:]:
    data.append(string.split())
data=np.asarray(data).astype(float)

x,y=data[:,0],data[:,1]

fig, ax = plt.subplots(nrows=2,ncols=2)

titles=['Density','Velocity x',r'Velocity $v_y$','Pressure P']
for i in range(2,6):

	z=data[:,i]
	
	'''
	RR = scipy.interpolate.RBFInterpolator(np.array([x, y]).T, z)
	zi = RR(np.array([x, y]).T)
	'''
	ax_actual=ax[int((i-2)/2),(i-2)%2]

	ax_actual.set_title(titles[i-2])
	ax_actual.axis('equal')
	#ax_actual.set_aspect('equal')
	#ax_actual.tricontour(x,y,zi,levels=512) 
	if (i==2):
		c=ax_actual.tricontourf(x,y,z,np.linspace(1, 3, 32),
                   	extend='both')
		cbar0=fig.colorbar(c,ax=ax_actual, orientation='vertical')
		cbar0.set_ticks(np.arange(1, 3.25, .25))
	
	elif (i==4):
		c=ax_actual.tricontourf(x,y,z,np.linspace(-0.12, 0.12, 32,endpoint=True),
                   	extend='both')
		cbar0=fig.colorbar(c,ax=ax_actual, orientation='vertical')
		cbar0.set_ticks(np.linspace(-0.12, 0.12, 6,endpoint=True))
	else:
		c=ax_actual.tricontourf(x,y,z,levels=32)
		cbar0=fig.colorbar(c,ax=ax_actual, orientation='vertical')

	
	'''
	z=z-z.min();
	z=z/z.max();
	z=z*y.max();
	ax_actual.plot(x,z,'o',markersize=1)
	'''
	ax_actual.set_xlim(x.min(),x.max())
	
	#ax_actual.set_ylim(0.5*y.max(),y.max())
	ax_actual.axis('equal')
	fig.tight_layout()


#name=sys.argv[2]
#plt.savefig("../../"+name+".png")
#plt.savefig("../../"+name+".svg")

#plt.gca().set_aspect('equal')
fig.tight_layout()	
plt.show()
"""
plt.imshow(zi, vmin=z.min(), vmax=z.max(), origin='lower',
           extent=[x.min(), x.max(), y.min(), y.max()])
plt.scatter(x, y, c=z)
plt.colorbar()
plt.show()
"""

