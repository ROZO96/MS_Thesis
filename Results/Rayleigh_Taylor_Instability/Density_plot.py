import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import sys


file_name=sys.argv[1]
file_name_split=file_name.split("/")
fp=open(file_name)
lines=fp.readlines()

n_point,time=int(lines[0].split()[0]),float(lines[0].split()[1])
print(n_point,time)



data = []
for string in lines[1:]:
    data.append(string.split())
data=np.asarray(data).astype(float)

x,y=data[:,0],data[:,1]
rho=data[:,2]
if 'vel' in file_name_split[1].split("_"):
	y=-y+2

ind_2=np.where((y>0.5*y.max()));

x=np.asarray([x[i] for i in ind_2]).flatten()
y=np.asarray([y[i] for i in ind_2]).flatten()
rho=np.asarray([rho[i] for i in ind_2]).flatten()


fig,ax=plt.subplots(ncols=2)
fig.suptitle('time=%1.3f'%(time))
c=ax[0].tricontour(x,y,rho,levels=[(rho.min()+rho.max())/2.0],colors='lightseagreen') 
#ax[0].clabel(c, inline=True, fontsize=10,colors='k')
ax[0].set_xlim(x.min(),x.max())
ax[0].set_ylim(y.min(),y.max())
ax[0].xaxis.set_major_locator(ticker.NullLocator())
ax[0].yaxis.set_major_locator(ticker.NullLocator())
ax[0].set_aspect('equal')

c=ax[1].tricontourf(x,y,rho,levels=32) 
c=ax[1].tricontourf(x,y,rho,levels=32) 
plt.colorbar(c)
ax[1].set_xlim(x.min(),x.max())
ax[1].set_ylim(y.min(),y.max())
ax[1].xaxis.set_major_locator(ticker.NullLocator())
ax[1].yaxis.set_major_locator(ticker.NullLocator())
ax[1].set_aspect('equal')

fig_name=file_name_split[0]+"/"+file_name_split[1]+"_density_plot_time_%.3f"%(time)+".svg"
print(fig_name)
plt.savefig(fig_name)



