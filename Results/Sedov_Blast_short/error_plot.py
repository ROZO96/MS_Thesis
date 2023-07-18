import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve, curve_fit
from tabulate import tabulate


print('Constructing analitical Solution....',end=' ')
##--------------------Analitical Solution---------------------##

#------Properties Before Shock Wave----#
rho_1=1;
P_1=100;
gamma=5.0/3.0;
R_blast=0.25;
E=1000/(gamma-1)#/(2*R_blast);
v1=np.sqrt(gamma*P_1/rho_1)
#------Problem Dimensionality-----#
nu=2; 

def sedov_sol(V):
	beta_2=(1-gamma)/(2*(gamma-1) +nu);
	beta_1= ((nu+2)*gamma/(2+nu*(gamma-1)))*(((2*nu*(2-gamma))/(gamma*((nu+2)**2)))- beta_2);
	beta_3=nu/(2*(gamma-1)+nu)
	beta_4=beta_1*(nu+2)/(2-gamma)
	beta_5=2/(gamma-2)
	beta_6=gamma/(2*(gamma-1)+nu)
	beta_7=(2+nu*(gamma-1))*beta_1/(nu*(2-gamma));

	lam=(((nu+2)*(gamma+1)/4)*V)**(-2/(2+nu));
	lam*=(((gamma+1)/(gamma-1))*((((nu+2)*gamma*V)/(2))-1))**(-beta_2);
	lam*=((((nu+2)*(gamma+1))/((gamma+1)*(nu+2)-2*(2+nu*(gamma-1))))*(1-((2+nu*(gamma-1))*V/2)))**(-beta_1);

	f=(nu+2)*(gamma+1)*(V/4)*lam;

	g=(((gamma+1)/(gamma-1))*(((nu+2)*gamma*V/2)-1))**beta_3;
	g*=(((gamma+1)/(gamma-1))*(1-((nu+2)*V/2)))**(beta_5);
	g*=((((nu+2)*(gamma+1))/((nu+2)*(gamma+1)-2*(2+nu*(gamma-1))))*(1-((2+nu*(gamma-1))*V/2)))**(beta_4);

	h=((nu+2)*(gamma+1)*V/4)**(2*nu/(2+nu));
	h*=(((gamma+1)/(gamma-1))*(1-((nu+2)*V/2)))**(beta_5+1);
	h*=((((nu+2)*(gamma+1))/((nu+2)*(gamma+1)-2*(2+nu*(gamma-1))))*(1-((2+nu*(gamma-1))*V/2)))**(beta_4-2*beta_1);
	
	
	return lam,g,f,h


def solve_V(V,lam):
	lam_resp=sedov_sol(V)[0];	
	return lam-lam_resp
	

def shock_update(r,t):
	#t+=to;
	#------Properties post Shock-----#
	r2=((E/rho_1)**(1/(2+nu))) *(t**(2/(2+nu)));# Shock Position
	v2=(4/((nu+2)*(gamma+1)))*((E/rho_1)**(1/2)) *(1/(r2**(nu/2)))
	rho2=((gamma+1)/(gamma-1))*rho_1
	p2=((8*E)/(((nu+2)**2)*(gamma+1)))*1/(r2**nu);
	#p2=P_1-(-rho2*(v2**2)+(v1**2)*rho_1)

	rho=[];
	v=[];
	p=[];
	V0=V_min;
	for r_i in r:
		if r_i<=r2:
			lam=(r_i)/r2;
			V0=fsolve(solve_V,V0,args=lam)[0]
			lam,g,f,h=sedov_sol(V0);		
			rho.append(rho2*g);
			v.append(v2*f);
			p.append(p2*h);
		else:
			rho.append(rho_1);
			v.append(0.00000001);
			p.append(P_1);	
	return rho,v,p
	
def cilindiral_conversion(x,y):
	CENTRE_X = 0.50*2;
	CENTRE_Y = 0.50*2;

	return np.sqrt((x-CENTRE_X)**2 +(y-CENTRE_Y)**2);
	
V_min=2/((nu+2)*gamma)




	
##------------------Results reading and error calculation----------------##
print('Calculating error for each Scheme')

scheme=['N','LDA','B','Bx']
case=[3,2,1,4,5]
variables=np.array([2,3,5])
variable_name=np.array([r'$\rho$',r'$V_r$'])
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
			r=np.asarray([r[i] for i in ind]).flatten()
			
			ind_2=np.where(r>=0.4);
			r=np.asarray([r[i] for i in ind_2]).flatten()
			
			exact_sol=shock_update(r,t);
		
		data_u=[]
		fig=plt.figure()
		for k in range(np.size(variables)):
			if k!=1:
				data_u_intial=np.asarray([data[i,variables[k]] for i in ind]).flatten()
				data_u.append(np.asarray([data_u_intial[i] for i in ind_2]).flatten())
			else:
				ux=np.asarray([data[i,variables[k]] for i in ind]).flatten()
				uy=np.asarray([data[i,variables[k]+1] for i in ind]).flatten()
				ux=np.asarray([ux[i] for i in ind_2]).flatten()
				uy=np.asarray([uy[i] for i in ind_2]).flatten()
				
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
fig,ax=plt.subplots(figsize=(30,10),nrows=1,ncols=np.size(variable_name))
colors=['red','green','blue','orange']
order=np.zeros((np.size(scheme),np.size(variable_name)))
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
print('Figure saved as Error_plot.svg in current directory'); 
print('\n Table of results for order of error estimation')
table_2=tabulate(np.c_[scheme,order],np.r_[np.array('Scheme'),variable_name], tablefmt="fancy_grid")
print(table_2)

with open('error_estimation_results.txt', 'w') as f:
    f.write(table_2)
print('Table saved as error_estimation_results.txt in current directory'); 
