import numpy as np
import matplotlib.pyplot as plt
import sys
from scipy.optimize import fsolve

file_directory=sys.argv[1]
print('Results for Sod Schock Tube in x Direction') 
gamma =5.0/3.0;
Gamma=(gamma-1)/(gamma+1)
beta=(gamma-1)/(2*gamma)
t=0.1



#------------------Region 4-------------------#
P4=1;
rho_4=1;
u4=0;
a4=np.sqrt(gamma*P4/rho_4);

#------------------Region 1-------------------#
P1=0.1;
rho_1=0.125
u1=0;
a1=np.sqrt(gamma*P1/rho_1);

def P_func(p3):
	u_2=(p3-P1)*np.sqrt((1-Gamma)/(rho_1*(p3+Gamma*P1))) 
	u_3=(P4**beta -p3**beta)*np.sqrt((1-Gamma**2)*(P4**(1/gamma))/(rho_4*Gamma**2))
	return u_2-u_3 

#------------------Region 2------------------#
P2=fsolve(P_func,P1)[0]
rho_2=rho_1*(P2+Gamma*P1)/(P2*Gamma+P1)
u2=(P2-P1)/np.sqrt((rho_1/2) *((gamma+1)*P2+(gamma-1)*P1))
#------------------Region 3------------------#
P3=P2;
rho_3=rho_4*(P3/P4)**(1/gamma)
u3=u2;

W=(a1*np.sqrt(((gamma+1)/(2*gamma))*((P3/P1) -1) +1) )

#------------------Intermediate Region 4-3------------------#
def properties_int(x):
	u=2/(gamma+1) *(a4+(x-2*0.75)/t)
	rho=rho_4*(1-(gamma-1)/2 * u/a4)** (2/(gamma-1))
	P=P4*(1-(gamma-1)/2 * u/a4)** (2*gamma/(gamma-1))
	return u,rho,P


def sod_profile(x_values):
	u=[]
	rho=[]
	P=[]
	for x in x_values:
		
			if x >=W*t+2*0.75:
				u.append(u1);
				rho.append(rho_1);
				P.append(P1)
			elif x>=u2*t+2*0.75:
				u.append(u2);
				rho.append(rho_2);
				P.append(P2);
			elif x>=(u3-np.sqrt(gamma*P3/rho_3))*t+2*0.75:
				u.append(u3);
				rho.append(rho_3);
				P.append(P3);
			elif x>=-a4*t+2*0.75:
				sol=np.array(properties_int(x))
				u.append(sol[0]);
				rho.append(sol[1]);
				P.append(sol[2]);
				
			else:
				u.append(u4);
				rho.append(rho_4);
				P.append(P4);
					
	return rho,u,P

print('Reading result files....')
file_B = open(file_directory+'/B_Scheme/snapshot_20.txt')
file_N = open(file_directory+'/N_Scheme/snapshot_20.txt')
file_LDA = open(file_directory+'/LDA_Scheme/snapshot_20.txt')
file_Bx = open(file_directory+'/Bx_Scheme/snapshot_20.txt')


lines_B=file_B.readlines();
lines_N=file_N.readlines();
lines_LDA=file_LDA.readlines();
lines_Bx=file_Bx.readlines();


n_point,t=int(lines_N[0].split()[0]),float(lines_N[0].split()[1])

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
ind_2=np.where((y<1.05) & (y>0.95));

x=np.asarray([x[i] for i in ind_2]).flatten()
y=np.asarray([y[i] for i in ind_2]).flatten()

ind = np.lexsort((y,x))

x=np.asarray([x[i] for i in ind])
y=np.asarray([y[i] for i in ind])

print('Plotting Results...')
variables=np.array([2,3,5])
variable_name=np.array([r'$\rho$',r'$U_x$',r'$P$'])
schemes=np.array(['N','LDA','B','Bx'])
exact_sol=sod_profile(x)

plt.rcParams.update({'font.size': 20})
fig,ax=plt.subplots(figsize=(30,10),nrows=1,ncols=3)

for j in range(np.size(variables)):
	ax[j].plot(x,exact_sol[j],label='Exact Sol.')
	ax[j].set_ylabel(variable_name[j])
	ax[j].set_xlabel(r'Position $x$')
	U=[]
	for scheme in schemes:
		if scheme=='N':
			U=np.asarray([data_N[i,variables[j]] for i in ind_2]).flatten()
		elif scheme=='LDA': 
			U=np.asarray([data_LDA[i,variables[j]] for i in ind_2]).flatten()
		elif scheme=='B':
			U=np.asarray([data_B[i,variables[j]] for i in ind_2]).flatten()
		elif scheme=='Bx':
			U=np.asarray([data_Bx[i,variables[j]] for i in ind_2]).flatten()
		
		U=[U[i] for i in ind]
		ax[j].plot(x,U,label=scheme)

	ax[j].set_xlim([1,2])
	ax[j].set_ylim([0,None])
	
ax[-1].legend(loc='upper right')
plt.tight_layout()
plt.savefig(file_directory+'/Results.svg')
print('Plot '+file_directory+ '/Results.svg saved.')
print('Done')

	
