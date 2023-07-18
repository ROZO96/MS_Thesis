import numpy as np 
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import sys
import glob

file_directory=sys.argv[1] 
files_B=glob.glob(file_directory+"/B_Scheme/sc_file_*")
files_Bx=glob.glob(file_directory+"/Bx_Scheme/sc_file_*")
bins=np.linspace(0,1,11)
labels=['B','Bx']
fig,ax=plt.subplots()
def animate(i):
	lines_B=open(files_B[i]).readlines();
	lines_Bx=open(files_Bx[i]).readlines();
	
	theta_B_values=[]
	theta_Bx_values=[]
	n_triangle,time=int(lines_B[0].split()[0]),float(lines_B[0].split()[1])
	for j in range(np.size(lines_B)-1):
		Theta_B = float(lines_B[j+1].split()[-1])	
		Theta_Bx = float(lines_Bx[j+1].split()[-1])
		theta_B_values.append(Theta_B);
		theta_Bx_values.append(Theta_Bx);
	weight_B = np.ones_like(theta_B_values)/float(len(theta_B_values))
	weight_Bx = np.ones_like(theta_Bx_values)/float(len(theta_Bx_values))
	results=np.c_[np.array(theta_B_values),np.array(theta_Bx_values)]
	weights=np.c_[weight_B ,weight_Bx]
	plt.clf() 
	plt.hist(results,bins=bins,label=labels,weights=weights)
	plt.title("time=%f"% (time));
	plt.ylim(0,1);
	plt.xlabel(r'$\Theta$')
	plt.ylabel(r'Frecuency Fraction')
	plt.legend();
	plt.draw();

ani = animation.FuncAnimation(fig, animate,save_count=np.size(files_B))
writergif = animation.PillowWriter(fps=30) 
ani.save(file_directory+'/Theta_histogram_evolution.gif', writer='pillow')

