import numpy as np 
import matplotlib.pyplot as plt
import sys

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

file_directory=sys.argv[1]
print('Results for 2D YEE Vortex')
print('Reading result files....')
file_B = open(file_directory+'/B_Scheme/snapshot_20.txt')
file_N = open(file_directory+'/N_Scheme/snapshot_20.txt')
file_LDA = open(file_directory+'/LDA_Scheme/snapshot_20.txt')
file_Bx = open(file_directory+'/Bx_Scheme/snapshot_20.txt')
'''
file_Bx_1 = open(file_directory+'/Bx_Scheme_1/snapshot_20.txt')
'''

lines_B=file_B.readlines();
lines_N=file_N.readlines();
lines_LDA=file_LDA.readlines();
lines_Bx=file_Bx.readlines();
'''
lines_Bx_1=file_Bx_1.readlines();
'''

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

'''
data_Bx_1 = []
for string in lines_Bx_1[1:]:
    data_Bx_1.append(string.split())
data_Bx_1=np.asarray(data_Bx_1).astype(float)
'''

x,y=data_N[:,0],data_N[:,1]
r=cilindiral_conversion(x,y);

ind = np.argsort(r)
r=np.asarray([r[i] for i in ind])
x=np.asarray([x[i] for i in ind])
y=np.asarray([y[i] for i in ind])
print('Plotting Results...');
variables=np.array([2,3,5]);
variable_name=np.array([r'$\rho$',r'$V_\theta$',r'$P$']);
schemes=np.array(['N','LDA','B','Bx']);
exact_sol=vortex_profile(r,x,y);

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
			'''
			elif scheme=='Bx_1':
				U=np.asarray([data_Bx_1[i,variables[j]] for i in ind]).flatten()
			'''
			
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
			'''
			elif scheme=='Bx_1':
				Ux=np.asarray([data_Bx_1[i,variables[j]] for i in ind]).flatten()
				Uy=np.asarray([data_Bx_1[i,variables[j]+1] for i in ind]).flatten()
				U=np.sqrt(Ux**2+Uy**2)
			'''
				
			ax[j].plot(r,U,label=scheme)
	ax[j].plot(r,exact_sol[j],label='Exact Sol.',color='k');
	ax[j].set_xlim([0,np.max(r)])
	#ax[j].set_ylim([0,None])
	
ax[-1].legend(loc='lower right')
plt.tight_layout()
plt.savefig(file_directory+'/Results.svg')
print('Plot '+file_directory+ '/Results.svg saved.')
print('Done')


