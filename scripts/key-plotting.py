#demonstration of plotfunctions. called from main UCLCHEM directory
#it reads full UCLCHEM output and saves a plot of the abudances of select species 

from plotfunctions import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.lines as mlines

#pick species, any number is fine
speciesNames=["0.25","0.5","1","2","4"]
input_file="output/full.dat"
plot_file="../VaryFromSolar/KeyAndOtherPlots1/KeyPlotbb.png"
varyfactor = [0.25, 0.5, 1, 2, 4]
#Linestles = [(0,(3,10,1,10)),(0,(3,5,1,5,1,5)),(0,()),(0,(3,5,1,5)),(0,(3,1,1,1))]
#Linestles = [(0,(4,2,1,2,1,2,1,2,1,2)),(0,(4,2,1,2,1,2)),(0,()),(0,(1,2,4,2,4,2)),(0,(1,2,4,2,4,2,4,2,4,2))]
Linestles = [(0,(1,2)),(0,(1,1)),(0,()),(0,(2,1,1,1)),(0,(2,1))]


#call read_uclchem. 
#time,dens,temp,abundances=read_uclchem(input_file,speciesNames)
time = [1.0,1e7]
abundances=[[1e-14,1e-14],[1e-12,1e-12],[1e-10,1e-10],[1e-8,1e-8],[1e-6,1e-6]]

#write out to columnated output,
#write_cols("output/democolumns.dat",time,dens,abundances)

#plot species and save to test.png, alternatively send dens instead of time.
fig,axs=plt.subplots(figsize=(16,9)) 
plt.subplots_adjust(right=0.8,left=None,bottom=None,top=None)

#plt.subplot(121)
rtists = []
factorlines = [mlines.Line2D([],[],color='w',label="factor",linestyle='None')]
cmapit = iter(cm.tab10(np.linspace(0,1,7)))
next(cmapit)
#plot_species will make an axis if you don't provide one
#it will save if plotFile is given otherwise just returns figure and axis

for j,v in enumerate(speciesNames):
    if j == 2:
        sn = "X"
    else:
        sn = "_X"
#    ax=plot_species(speciesNames[j],time,abundances[j],ax=axs,ls=Linestles[j],lw=3.0)
    rtist = axs.plot(time,abundances[j],color="r",label=sn,linestyle=Linestles[j],linewidth=3.0)
    rtists.append(rtist)
    factorlines.append( mlines.Line2D([],[],color='k',label=v,linestyle=Linestles[j]))
#ax2=ax.twinx()
#ax2.plot(time,temp)
lbls = ["species"]
hdls = [mlines.Line2D([],[],color='w',label="species",linestyle='None')]
hndls,labs = axs.get_legend_handles_labels()
lbls.extend(labs)
lbls.append("factor")
lbls.extend(speciesNames)
hdls.extend(hndls)
hdls.extend(factorlines)


#plot species returns the axis so we can further edit
axs.set(xscale='linear',yscale='log',xlim=(1,1e7),ylim=(5e-18,3e-3))
axs.set_ylabel('X/H',fontsize='x-large')
axs.set_xlabel('Time / Myears',fontsize='x-large')
axs.ticklabel_format(axis='x',useMathText=True)
axs.set_title('Key for 5 levels of varied element',fontsize='xx-large')


axs.legend(hdls,lbls,loc='upper left',bbox_to_anchor=(1.0,1.0),fontsize='xx-large',handlelength=5,borderaxespad=0.0)#loc='upper left',bbox_to_anchor=(0.0,0.0,1.0,1.0),
#plt.show()
#overwrite our previous plot
fig.savefig(plot_file)
