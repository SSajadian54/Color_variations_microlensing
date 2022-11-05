import numpy as np
import matplotlib.pyplot as plt
import pylab as py
from matplotlib import rcParams
from numpy import ma
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
rcParams["font.size"] = 11.5
rcParams["font.family"] = "sans-serif"
rcParams["font.sans-serif"] = ["Computer Modern Sans"]
rcParams["text.usetex"] = True
rcParams["text.latex.preamble"] = r"\usepackage{cmbright}"
import math
from mpl_toolkits.axes_grid.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
import matplotlib as mpl
from matplotlib import gridspec
from matplotlib.ticker import FormatStrFormatter
#################################################################################3
n=[0,0,0,0]
n[1]=28348
n[0]=n[2]=31494
sou=np.zeros((4,n[0],4))
sou[0,:,:]=np.loadtxt("LimbD.dat")
sou[1,:n[1],:]=np.loadtxt("LimbG.dat")
sou[2,: ,:]=np.loadtxt("Spot.dat")
f1=open("./sou_4_2_90.dat","r")
ns= sum(1 for line in f1)
sou1=np.zeros((ns,12))
sou1= np.loadtxt("./sou_4_2_90.dat")
k=int(0)
for h in range(ns):  
    if(sou1[h,1]>0.0):
        sou[3,k,0]=sou1[h,2]*(0.99/1.88609)    
        sou[3,k,1]=sou1[h,3]*(0.99/1.88609)      
        sou[3,k,2]=sou1[h,6]        
        k+=1
n[3]=k      
print n[0], n[1], n[2], n[3]

        

fig, axs = plt.subplots(nrows=2, ncols=2,figsize=(8,8))
m=0; 
for ax in axs.flat:
    fig.subplots_adjust(wspace=0,hspace=0)
    b= int(n[m])
    ax.scatter(sou[m,:b,0],sou[m,:b,1],c=sou[m,:b,2],cmap='viridis',marker='o',s=25,linewidth=0)  
    ax.set_xlim(-1.1,1.1)
    ax.set_ylim(-1.1,1.1)
    ax.set_aspect('equal', adjustable='box')            
    py.setp(ax.get_xticklabels(), visible=False)
    py.setp(ax.get_yticklabels(), visible=False)
    plt.gca().axes.get_yaxis().set_visible(False)
    ax.axis('off')
    if(m==0):  ax.set_title(r"$\rm{Limb-Darkening}$",fontsize=16)
    if(m==1):  ax.set_title(r"$\rm{Gravity-Darkening}$",fontsize=16)
    if(m==2):  ax.text(-.5,-1.2,r"$\rm{Stellat~Spot}$",fontsize=16,rotation=0)
    if(m==3):  ax.text(-.7,-1.2,r"$\rm{Non-Radial~Pulsating}$",fontsize=16,rotation=0)
    m+=1
all_axes = fig.get_axes()
plt.savefig("source.jpg",dpi=200)

    
            
                


'''

plt.clf()
fig= plt.figure()
plt.scatter(sou2[:,0],sou2[:,1],c=sou2[:,2],cmap='viridis', marker='o',s=30.0,linewidth=0)      
py.xlim([-1.1,1.1])
py.ylim([-1.1,1.1])
plt.gca().set_aspect('equal', adjustable='box')
plt.savefig("limbD.jpg",dpi=200)

plt.clf()
fig= plt.figure()
plt.scatter(sou2[:,0],sou2[:,1],c=sou2[:,2],cmap='viridis', marker='o',s=30.0,linewidth=0)      
py.xlim([-1.1,1.1])
py.ylim([-1.1,1.1])
plt.gca().set_aspect('equal', adjustable='box')
plt.savefig("limbG.jpg",dpi=200)

########################################################################

sou2=np.zeros((31494,4))
sou2=np.loadtxt("Spot.dat")
plt.clf()
fig= plt.figure()
plt.scatter(sou2[:,0],sou2[:,1],c=sou2[:,2],cmap='viridis', marker='o',s=30.0,linewidth=0)      
py.xlim([-1.1,1.1])
py.ylim([-1.1,1.1])
plt.gca().set_aspect('equal', adjustable='box')
plt.savefig("Spot.jpg",dpi=200)


########################################################################

plt.clf()
fig= plt.figure()
plt.scatter(sou2[:k2,2],sou2[:k2,3],c=sou2[:k2,6],cmap='viridis', marker='o',s=25.0,linewidth=0)      
py.xlim([-2.0,2.0])
py.ylim([-2.0,2.0])
plt.gca().set_aspect('equal', adjustable='box')
plt.savefig("puls.jpg",dpi=200)

'''  
        
