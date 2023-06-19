import numpy as np
import matplotlib.pyplot as plt


file_B = open('B Scheme/snapshot_20.txt')
file_N = open('N Scheme/snapshot_20.txt')
file_LDA = open('LDA Scheme/snapshot_20.txt')
file_Bx = open('Bx Scheme/snapshot_20.txt')

lines_B=file_B.readlines();
lines_N=file_N.readlines();
lines_LDA=file_LDA.readlines();
lines_Bx=file_Bx.readlines();


n_point,time=int(lines_N[0].split()[0]),float(lines_N[0].split()[1])
#print(n_point,time)

data_N = []
for string in lines_N[1:]:
    data_N.append(string.split())
data_N=np.asarray(data_N).astype(float)


data_LDA = []
for string in lines_LDA[1:]:
    data_LDA.append(string.split())
data_LDA=np.asarray(data_LDA).astype(float)

data_B = []
for string in lines_B[1:]:
    data_B.append(string.split())
data_B=np.asarray(data_B).astype(float)

data_Bx = []
for string in lines_Bx[1:]:
    data_Bx.append(string.split())
data_Bx=np.asarray(data_Bx).astype(float)

x,y=data_N[:,0],data_N[:,1]
ind_2=np.where((y<1.025) & (y>0.975));

x=np.asarray([x[i] for i in ind_2]).flatten()
y=np.asarray([y[i] for i in ind_2]).flatten()

rho_N=np.asarray([data_N[i,2] for i in ind_2]).flatten()
rho_LDA=np.asarray([data_LDA[i,2] for i in ind_2]).flatten()
rho_B=np.asarray([data_B[i,2] for i in ind_2]).flatten()
rho_Bx=np.asarray([data_Bx[i,2] for i in ind_2]).flatten()

ind = np.lexsort((y,x))

x=np.asarray([x[i] for i in ind])
y=np.asarray([y[i] for i in ind])
rho_N=[rho_N[i] for i in ind]
rho_LDA=[rho_LDA[i] for i in ind]
rho_B=[rho_B[i] for i in ind]
rho_Bx=[rho_Bx[i] for i in ind]


xi, yi = np.meshgrid(x,y)
"""
fig, ax = plt.subplots(nrows=2,ncols=2)
ax_N=ax[0,0];
ax_LDA=ax[0,1];
ax_B=ax[1,0];
ax_Bx=ax[1,1];

ax_N.set_title('N')
ax_LDA.set_title('LDA')
ax_B.set_title('B')
ax_Bx.set_title('Bx')

ax_N.plot(x,rho_N)
ax_LDA.plot(x,rho_LDA)
ax_B.plot(x,rho_B)
ax_Bx.plot(x,rho_Bx)

ax_N.set_xlim(x.min(),x.max())
ax_LDA.set_xlim(x.min(),x.max())
ax_B.set_xlim(x.min(),x.max())
ax_Bx.set_xlim(x.min(),x.max())


fig.tight_layout()
"""
plt.figure()

plt.plot(x,rho_N,label='N')
plt.plot(x,rho_LDA,label='LDA')
plt.plot(x,rho_B,label='B')
plt.plot(x,rho_Bx,label='Bx')
plt.xlabel(r'Position $x$')
plt.ylabel(r'Density $\rho$')
plt.legend()

plt.show()






