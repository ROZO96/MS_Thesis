import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve, curve_fit
from tabulate import tabulate
import warnings
warnings.filterwarnings("ignore")

print('Constructing analitical Solution....',end=' ')
##--------------------Analitical Solution---------------------##
def P_func(p3):
	u_2=(p3-P1)*np.sqrt((1-Gamma)/(rho_1*(p3+Gamma*P1))) 
	u_3=(P4**beta -p3**beta)*np.sqrt((1-Gamma**2)*(P4**(1/gamma))/(rho_4*Gamma**2))
	return u_2-u_3 

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
#-------------Functions-----------------#
gamma =5.0/3.0;
Gamma=(gamma-1)/(gamma+1)
beta=(gamma-1)/(2*gamma)

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

print('Done \n') 

##------------------Results reading and error calculation----------------##
print('Calculating error for each Scheme')

scheme=['N','LDA','B','Bx']
case=[3,2,1,4,5]
variables=np.array([2,3,5])
variable_name=np.array([r'$\rho$',r'$U_x$',r'$P$'])
error_data=np.zeros((np.size(scheme),np.size(case),np.size(variables)))
h=np.zeros(np.size(case))
for l in range(np.size(case)):
	x,y=[],[];
	print('\n')
	print('Case', case[l])
	ind=[];
	ind_2=[];
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
			
			ind_2=np.where(x>=1);
			x=np.asarray([x[i] for i in ind_2]).flatten()
			y=np.asarray([y[i] for i in ind_2]).flatten()
			
			ind_3=np.where((y<2.05) & (y>-0.95));
			
			x=np.asarray([x[i] for i in ind_3]).flatten()
			y=np.asarray([y[i] for i in ind_3]).flatten()

			ind = np.lexsort((y,x))

			x=np.asarray([x[i] for i in ind])
			y=np.asarray([y[i] for i in ind])
			
			exact_sol=sod_profile(x)
		
		data_u=[]
		fig=plt.figure()
		for k in range(np.size(variables)):
			data_u_intial=np.asarray([data[i,variables[k]] for i in ind_2]).flatten()
			data_u_intial=np.asarray([data_u_intial[i] for i in ind_3]).flatten()
			data_u.append(np.asarray([data_u_intial[i] for i in ind]).flatten())
		data_u=np.array(data_u)
		#error=np.max(np.abs(data_u-exact_sol),axis=1)
		error=np.sum(np.abs(data_u-exact_sol),axis=1)/np.shape(data_u)[1]
		#error=np.sqrt(np.sum((data_u-exact_sol)**2,axis=1))/np.shape(data_u)[1]
		
		
		error_data[j,l]=error;
		print('Done')
##------------------Error ploting and order of the error calculation----------------##
def powlaw(x, a, b) :
    return a * (x**b)
print('\n Creating Error Plot.....',end=' ')
plt.rcParams.update({'font.size': 20})
fig,ax=plt.subplots(figsize=(30,10),nrows=1,ncols=3)
colors=['red','green','blue','orange']
order=np.zeros((np.size(scheme),np.size(variables)))
for i in range(np.size(scheme)):
	for j in range(np.size(variable_name)):
		popt, pcov=curve_fit(powlaw, h, error_data[i,:,j], maxfev=2000)
		z = np.polyfit(np.log(h), np.log(error_data[i,:,j]), 1)
		print(z)
		p = np.poly1d(z)

		order[i,j]=popt[1]
		ax[j].set_ylabel('Error '+ variable_name[j])
		ax[j].set_xlabel(r'$h$')
		ax[j].loglog(h,error_data[i,:,j],'.-',color=colors[i])
		#ax[j].loglog(h,powlaw(h,*popt),'--',color=colors[i],label=scheme[i])	
		ax[j].loglog(h,np.exp(p(np.log(h))),'--',color=colors[i],label=scheme[i])
ax[-1].legend(loc='lower right')		
plt.tight_layout()
print('Done')
plt.savefig('Error_plot.svg')
print('Figure saved as Error_plot.svg in current directory'); 
print('\n Table of results for order of error estimation')
table_2=tabulate(np.c_[scheme,order],np.r_[np.array('Scheme'),variable_name], tablefmt="fancy_grid")
print(table_2)

with open('error_estimation_results.txt', 'w') as f:
    f.write(table_2)
print('Table saved as error_estimation_results.txt in current directory'); 
		
		
		
		
