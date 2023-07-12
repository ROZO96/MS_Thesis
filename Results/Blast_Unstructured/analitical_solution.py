import numpy as np 
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import sys


#------Properties Before Shock Wave----#
rho_1=1;
P_1=100;
gamma=5.0/3.0;
R_blast=0.25;
E=1000#/(2*R_blast);
v1=np.sqrt(gamma*P_1/rho_1)
#------Problem Dimensionality-----#
nu=2; 

'''
def solve_t_initial(t):
	r2=((E/rho_1)**(1/(2+nu))) *(t**(2/(2+nu)));
	return r2-R_blast;

to=fsolve(solve_t_initial,1e-6)[0]
'''


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
	#p2=((8*E)/(((nu+2)**2)*(gamma+1)))*1/(r2**nu);
	p2=P_1-(-rho2*(v2**2)+(v1**2)*rho_1)

	rho=[];
	v=[];
	p=[];
	V0=V_min;
	for r_i in r:
		if (r_i-R_blast)<0:
			lam,g,f,h=sedov_sol(V_min);		
			rho.append(rho2*g);
			v.append(v2*f);
			p.append(p2*h);
		elif r_i-R_blast<=r2:
			lam=(r_i-R_blast)/r2;
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
	CENTRE_X = 0.50*10;
	CENTRE_Y = 0.50*10;

	return np.sqrt((x-CENTRE_X)**2 +(y-CENTRE_Y)**2);
	
V_min=2/((nu+2)*gamma)


file_directory=sys.argv[1]
print('Results for 2D Sedov Blast')
print('Reading result files....')
file_B = open(file_directory+'/B_Scheme/snapshot_20.txt')
file_N = open(file_directory+'/N_Scheme/snapshot_20.txt')
file_LDA = open(file_directory+'/LDA_Scheme/snapshot_20.txt')
file_Bx = open(file_directory+'/Bx_Scheme/snapshot_20.txt')

lines_B=file_B.readlines();
lines_N=file_N.readlines();
lines_LDA=file_LDA.readlines();
lines_Bx=file_Bx.readlines();

n_point,time=int(lines_N[0].split()[0]),float(lines_N[0].split()[1])	


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
r=cilindiral_conversion(x,y);


ind = np.argsort(r)
r=np.asarray([r[i] for i in ind])
print('Plotting Results...');
variables=np.array([2,3,5]);
variable_name=np.array([r'$\rho$',r'$V_r$',r'$P$']);
schemes=np.array(['N','LDA','B','Bx']);
exact_sol=shock_update(r,time);


plt.rcParams.update({'font.size': 20})
fig,ax=plt.subplots(figsize=(30,10),nrows=1,ncols=3)
for j in range(np.size(variables)):
	
	ax[j].set_ylabel(variable_name[j]);
	ax[j].set_xlabel(r'Position $r$');
	U=[];
	for scheme in schemes:
		if j!=1: 
			if scheme=='N':
				U=np.asarray([data_N[i,variables[j]] for i in ind]).flatten()
			elif scheme=='LDA': 
				U=np.asarray([data_LDA[i,variables[j]] for i in ind]).flatten()
			elif scheme=='B':
				U=np.asarray([data_B[i,variables[j]] for i in ind]).flatten()
			elif scheme=='Bx':
				U=np.asarray([data_Bx[i,variables[j]] for i in ind]).flatten()
		
			ax[j].plot(r,U,label=scheme)
		else:
			if scheme=='N':
				Ux=np.asarray([data_N[i,variables[j]] for i in ind]).flatten()
				Uy=np.asarray([data_N[i,variables[j]+1] for i in ind]).flatten()
				U=np.sqrt(Ux**2+Uy**2);
			elif scheme=='LDA': 
				Ux=np.asarray([data_LDA[i,variables[j]] for i in ind]).flatten()
				Uy=np.asarray([data_LDA[i,variables[j]+1] for i in ind]).flatten()
				U=np.sqrt(Ux**2+Uy**2);
			elif scheme=='B':
				Ux=np.asarray([data_B[i,variables[j]] for i in ind]).flatten()
				Uy=np.asarray([data_B[i,variables[j]+1] for i in ind]).flatten()
				U=np.sqrt(Ux**2+Uy**2)
			elif scheme=='Bx':
				Ux=np.asarray([data_Bx[i,variables[j]] for i in ind]).flatten()
				Uy=np.asarray([data_Bx[i,variables[j]+1] for i in ind]).flatten()
				U=np.sqrt(Ux**2+Uy**2)
				
			ax[j].plot(r,U,label=scheme)
	ax[j].plot(r,exact_sol[j],label='Exact Sol.',color='k');
	ax[j].set_xlim([0,np.max(r)])
	#ax[j].set_ylim([0,None])
	
ax[-1].legend(loc='lower right')
plt.tight_layout()
plt.savefig(file_directory+'/Results.svg')
print('Plot '+file_directory+ '/Results.svg saved.')
print('Done')

	
	
