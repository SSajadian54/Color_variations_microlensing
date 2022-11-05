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
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)
from matplotlib.ticker import FixedLocator, FixedFormatter
#######################################################################################

f0= open("./param_all.txt","r")
nn= sum(1 for line in f0)
print "No.  row in param1.dat file:  ",  nn
par= np.zeros((int(nn),26))
par=np.loadtxt("param_all.txt") 

nam=[r"$u_{0}(\rho_{\star})$",r"$,~\rm{t_{E}(day)}$",r"$,~\rm{T_{\rm{eff}}(K)}$", r"$,~\rm{\log_{10}[g]}$"]


for n in range(80):

    
    m=n
    u0= float(par[m,10]); tE=float(par[m,9]);  tstar= float(par[m,12])
    rho=float(par[m,22]); teff= int(par[m, 20]);  logg=float(par[m,18])
    print "**************** No. of light curve:  ", n
    print "  u0:  ", u0,  "   tE: ", tE, "rho:  ", rho, "  tef:  ", teff, "   logg:  ", logg
    f1=open("./MONT/M_{0:d}.dat".format(n),"r")
    nm= sum(1 for line in f1)
    print "No.  of mag.dat:  ",  nm
    plt.clf()
    plt.cla()
    if(nm>0):
        mag=np.zeros((nm,10)); 
        mag=np.loadtxt("./MONT/M_{0:d}.dat".format(n)) 
        ma=np.zeros((7,2))
        ma[0,0]=np.min(mag[:,3]); ma[0,1]= np.max(mag[:,3]);
        ma[1,0]=np.min(mag[:,4]); ma[1,1]= np.max(mag[:,4]);
        ma[2,0]=np.min(mag[:,5]); ma[2,1]= np.max(mag[:,5]);
        ma[3,0]=np.min(mag[:,6]); ma[3,1]= np.max(mag[:,6]);
        ma[4,0]=np.min(mag[:,7]); ma[4,1]= np.max(mag[:,7]);
        miny= np.min(ma[:5,0]);   maxy= np.max(ma[:5,1]);  
        dy= abs((maxy-miny)*0.05)       
        minx=-5.0;   maxx=5.0;  dx=1.0;
        print "miny: ", miny, "   maxy:  ", maxy
        ma[0,0]=np.min(mag[:,8]);  ma[0,1]= np.max(mag[:,8]);
        ma[1,0]=np.min(mag[:,9]);  ma[1,1]= np.max(mag[:,9]);
        #ma[2,0]=np.min(mag[:,16]); ma[2,1]= np.max(mag[:,16]);
        #ma[3,0]=np.min(mag[:,17]); ma[3,1]= np.max(mag[:,17]);
        #ma[4,0]=np.min(mag[:,18]); ma[4,1]= np.max(mag[:,18]);
        #ma[5,0]=np.min(mag[:,19]); ma[5,1]= np.max(mag[:,19]);
        miny2= np.min(ma[:2,0]);   maxy2= np.max(ma[:2,1]);  
        print "miny2: ", miny2, "   maxy2:  ", maxy2
        dy2= abs((maxy2-miny2)*0.05)  
        
        

    plt.clf()
    #fig=plt.figure()
    gs = gridspec.GridSpec(2,1, height_ratios=[2, 1])
    axs0=plt.subplot(gs[0])
    axs0.plot(mag[:,0], mag[:,3], c="m", lw=1.)
    axs0.plot(mag[:,0], mag[:,4], c="blue", lw=1.)
    axs0.plot(mag[:,0], mag[:,5], c="green", lw=1.)
    axs0.plot(mag[:,0], mag[:,6], c="orange", lw=1.)
    axs0.plot(mag[:,0], mag[:,7], c="red", lw=1.)
    axs0.set_yticks(np.arange(miny,maxy,dy*3))
    axs0.set_xticks(np.arange(minx,maxx,dx))
    axs0.xaxis.set_major_formatter(FormatStrFormatter('%d'))
    axs0.yaxis.set_major_formatter(FormatStrFormatter('%d'))
    axs0.set_ylim(miny,maxy)
    axs0.set_xlim(minx,maxx)
    axs0.set_title(str(nam[0])+ '={0:.2f}'.format(u0) + str(nam[2])+ '={0:d}'.format(teff) + str(nam[3])+ '={0:.1f}'.format(logg),fontsize=14,color='k')
    axs0.set_ylabel(r"$\rm{Magnification}$", fontsize=16.)
    plt.text(-4.5,maxy-2*dy,r"$\rm{U-band}$",color="m", fontsize=16)
    plt.text(-4.5,maxy-4*dy,r"$\rm{B-band}$",color="b", fontsize=16)
    plt.text(-4.5,maxy-6*dy,r"$\rm{V-band}$",color="green", fontsize=16)
    plt.text(-4.5,maxy-8*dy,r"$\rm{R-band}$",color="orange", fontsize=16)
    plt.text(-4.5,maxy-10*dy,r"$\rm{I-band}$",color="red", fontsize=16)
  
  
  
    axs1= plt.subplot(gs[1],sharex= axs0)
    axs1.plot(mag[:,0], mag[:,8],"k-", lw=1.)
    axs1.plot(mag[:,0], mag[:,9],"grey", lw=1.)
    axs1.set_yticks(np.arange(miny2,maxy2,dy2*5.0))
    axs1.set_ylim(miny2,maxy2)
    axs1.set_xlim(minx,maxx)
    axs1.set_xticks(np.arange(minx,maxx,dx))
    axs1.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
    axs1.xaxis.set_major_formatter(FormatStrFormatter('%d'))
    axs1.set_xlabel(r"$\rm{(t-t_{0})/t_{\star}}$", fontsize=17,labelpad=0.1)
    axs1.set_ylabel(r"$\rm{\Delta m(mag)}$",  fontsize=17,labelpad=0.1)
    plt.setp(axs0.get_xticklabels(), visible=False)
    yticks = axs1.yaxis.get_major_ticks()
    plt.text(-4.5,maxy2-3*dy2,r"$\rm{V-I}$",color="black", fontsize=16)
    plt.text(-4.5,maxy2-6*dy2,r"$\rm{B-R}$",color="grey", fontsize=16)
    #yticks[-1].label1.set_visible(False)
    mpl.rcParams['axes.formatter.use_mathtext'] = 'True'
    mpl.rcParams['axes.formatter.useoffset'] = 'False'
    mpl.rcParams['figure.subplot.left'] = 1.0
    mpl.rcParams['figure.subplot.right'] = .95
    mpl.rcParams['figure.subplot.bottom'] = .20
    mpl.rcParams['figure.subplot.top'] = .90
    plt.subplots_adjust(hspace=.0)
    fig=plt.gcf()
    fig.savefig("./lightcurves/LI_{0:d}.jpg".format(n),dpi=200)
    print "First light curve was plotted ***** "
    #input("Enter a number ")
    py.clf()
############################################################################



