import  numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from matplotlib.widgets import Slider, Button, RadioButtons

#------Properties Before Shock Wave----#
rho_1=1;
P_1=1e-8;
gamma=5.0/3.0;
R_blast=0.25;
E=10#/(np.pi*R_blast**2);

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
	#------Properties post Shock-----#
	r2=((E/rho_1)**(1/(2+nu))) *(t**(2/(2+nu)));# Shock Position
	v2=(4/((nu+2)*(gamma+1)))*((E/rho_1)**(1/2)) *(1/(r2**(nu/2)))
	rho2=((gamma+1)/(gamma-1))*rho_1
	p2=((8*E)/(((nu+2)**2)*(gamma+1)))*1/(r2**nu)
	
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
		elif (r_i-R_blast)<=r2:
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
	
def max_values(t):
#------Properties post Shock-----#
	r2=((E/rho_1)**(1/(2+nu))) *(t**(2/(2+nu)));# Shock Position
	v2=(4/((nu+2)*(gamma+1)))*((E/rho_1)**(1/2)) *(1/(r2**(nu/2)))
	rho2=((gamma+1)/(gamma-1))*rho_1
	p2=((8*E)/(((nu+2)**2)*(gamma+1)))*1/(r2**nu);
	
	return rho2,v2,p2
	
V_min=2/((nu+2)*gamma)
'''
V_max=4/((nu+2)*(gamma+1));	
V=np.linspace(V_min,V_max,10000)
lam,g,f,h=sedov_sol(V);'''

fig,ax=plt.subplots(figsize=(30,10),nrows=1,ncols=3)
rmax=5*np.sqrt(2)
r=np.linspace(0,rmax,1001);
time_0=1e-6;
time_1=10;
fig.subplots_adjust(bottom=0.25)

rho,v,p=shock_update(r,time_0);

ax[0].plot(r,rho)
ax[1].plot(r,v)
ax[2].plot(r,p)
rho_max,v_max,p_max=max_values(time_0)

ax[0].set_ylim([None,rho_max])
ax[1].set_ylim([None,v_max])
ax[2].set_ylim([None,p_max])

ax[-1].legend(loc='upper right')

axtime= fig.add_axes([0.25, 0.1, 0.65, 0.03])
time_slider = Slider(
    ax=axtime,
    label='Time',
    valmin=time_0,
    valmax=time_1,
    valinit=time_0,
)
# The function to be called anytime a slider's value changes
def update(val):
	rho,v,p=shock_update(r,val);
	ax[0].cla()
	ax[1].cla()
	ax[2].cla()
	ax[0].plot(r,rho)
	ax[1].plot(r,v)
	ax[2].plot(r,p)
	#ax[0].set_ylim([None,rho_max])
	#ax[1].set_ylim([None,v_max])
	#ax[2].set_ylim([None,p_max])
	fig.canvas.draw_idle()
	

time_slider.on_changed(update)

plt.show()





