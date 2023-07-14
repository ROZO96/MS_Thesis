import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve, curve_fit
from tabulate import tabulate


print('Constructing analitical Solution....',end=' ')
##--------------------Analitical Solution---------------------##
def vortex_profile(r_values,x_values,y_values):
	CENTRE_X = 0.50*10;
	CENTRE_Y = 0.50*10;
	RHO_INF=1.0;
	PRESSURE_INF=1.0;
	K=1.0;
	GAMMA=5.0/3.0;
	BETA=5.0;
	
	P_LIST=[];
	V_LIST=[];
	RHO_LIST=[];
	for R,X,Y in zip(r_values,x_values,y_values):
		OMEGA_R=(BETA/(2.0*np.pi))*np.exp((1.0-R*R)/2.0);
		X_VEL=-OMEGA_R*(Y-CENTRE_Y);
		Y_VEL=OMEGA_R*(X-CENTRE_X);
		T_R=(PRESSURE_INF/RHO_INF)-((GAMMA-1)/GAMMA)*(BETA*BETA/(8.0*np.pi*np.pi))*np.exp(1.0-R*R);
		
		RHO=(T_R/K)**(1/(GAMMA-1));
		PRESSURE=K*(RHO**GAMMA);
		
		V_THETA=np.sqrt((X_VEL**2)+(Y_VEL**2));
		
		P_LIST.append(PRESSURE);
		V_LIST.append(V_THETA);
		RHO_LIST.append(RHO);
		
	P_LIST=np.array(P_LIST);
	V_LIST=np.array(V_LIST);
	RHO_LIST=np.array(RHO_LIST);
	
	return RHO_LIST, V_LIST, P_LIST;

def cilindiral_conversion(x,y):
	CENTRE_X = 0.50*10;
	CENTRE_Y = 0.50*10;

	return np.sqrt((x-CENTRE_X)**2 +(y-CENTRE_Y)**2);
	
	
##------------------Results reading and error calculation----------------##
print('Calculating error for each Scheme')

scheme=['N','LDA','B','Bx']
case=[5,4,1,2,3]
variables=np.array([2,3,5])
variable_name=np.array([r'$\rho$',r'$V_\theta$',r'$P$'])
error_data=np.zeros((np.size(scheme),np.size(case),np.size(variables)))
h=np.zeros(np.size(case))
for l in range(np.size(case)):
	x,y=[],[];
	print('\n')
	print('Case', case[l])
	ind=[];
	for j in range(np.size(scheme)):
		print('Scheme',scheme[j])
		print('Reading result files....',end=' ')
		file_f= open('Case_'+str(case[l])+'/'+scheme[j]+'_Scheme/snapshot_20.txt')
		lines=file_f.readlines();
		n_point,t=int(lines[0].split()[0]),float(lines[0].split()[1]);
		data=[]
		for string in lines[1:]:
    			data.append(string.split())
		data=np.asarray(data).astype(float)
		
		if (j==0):
			x,y=data[:,0],data[:,1]
			h[l]=np.max(x)/np.sqrt(np.size(x));
			r=cilindiral_conversion(x,y);
			
			ind = np.argsort(r)
			r=np.asarray([r[i] for i in ind])
			x=np.asarray([x[i] for i in ind])
			y=np.asarray([y[i] for i in ind])
			exact_sol=vortex_profile(r,x,y)
		
		data_u=[]
		fig=plt.figure()
		for k in range(np.size(variables)):
			if k!=1:
				data_u.append(np.asarray([data[i,variables[k]] for i in ind]).flatten())
			else:
				ux=np.asarray([data[i,variables[k]] for i in ind]).flatten()
				uy=np.asarray([data[i,variables[k]+1] for i in ind]).flatten()
				data_u.append(np.sqrt(ux**2+uy**2));
		data_u=np.array(data_u)
		#error=np.max(np.abs(data_u-exact_sol),axis=1)
		error=np.sum(np.abs(data_u-exact_sol),axis=1)/np.shape(data_u)[1]
		#error=np.sqrt(np.sum((data_u-exact_sol)**2,axis=1))/np.shape(data_u)[1]
		
		
		error_data[j,l]=error;
		print('Done')

##------------------Error ploting and order of the error calculation----------------##
print('\n Creating Error Plot.....',end=' ')
plt.rcParams.update({'font.size': 20})
fig,ax=plt.subplots(figsize=(30,10),nrows=1,ncols=3)
colors=['red','green','blue','orange']
order=np.zeros((np.size(scheme),np.size(variables)))
for i in range(np.size(scheme)):
	for j in range(np.size(variable_name)):
		z = np.polyfit(np.log(h), np.log(error_data[i,:,j]), 1)
		p = np.poly1d(z);
		order[i,j]=z[0]
		ax[j].set_ylabel('Error '+ variable_name[j])
		ax[j].set_xlabel(r'$h$')
		ax[j].loglog(h,error_data[i,:,j],'.-',color=colors[i])
		ax[j].loglog(h,np.exp(p(np.log(h))),'--',color=colors[i],label=scheme[i])
ax[-1].legend(loc='lower right')		
plt.tight_layout()
print('Done')
plt.savefig('Error_plot.svg')
plt.close()
print('Figure saved as Error_plot.svg in current directory'); 
print('\n Table of results for order of error estimation')
table_2=tabulate(np.c_[scheme,order],np.r_[np.array('Scheme'),variable_name], tablefmt="fancy_grid")
print(table_2)

with open('error_estimation_results.txt', 'w') as f:
    f.write(table_2)
print('Table saved as error_estimation_results.txt in current directory'); 
