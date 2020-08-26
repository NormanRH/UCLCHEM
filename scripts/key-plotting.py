#demonstration of plotfunctions. called from main UCLCHEM directory
#it reads full UCLCHEM output and saves a plot of the abudances of select species 

from plotfunctions import *

#pick species, any number is fine
speciesNames=[["X/4"],["X/2"],["Xx1"],["Xx2"],["Xx4"]]
input_file="output/full.dat"
plot_file="../VaryFromSolar/KeyAndOtherPlots1/KeyPlot.png"
varyfactor = [0.25, 0.5, 1, 2, 4]
Linestles = [(0,(3,10,1,10)),(0,(3,5,1,5,1,5)),(0,()),(0,(3,5,1,5)),(0,(3,1,1,1))]


#call read_uclchem. 
#time,dens,temp,abundances=read_uclchem(input_file,speciesNames)
time = [1.0,1e7]
abundances=[[[1e-14,1e-14]],[[1e-12,1e-12]],[[1e-10,1e-10]],[[1e-8,1e-8]],[[1e-6,1e-6]]]

#write out to columnated output,
#write_cols("output/democolumns.dat",time,dens,abundances)

#plot species and save to test.png, alternatively send dens instead of time.
fig,axs=plt.subplots(figsize=(16,9)) 

#plot_species will make an axis if you don't provide one
#it will save if plotFile is given otherwise just returns figure and axis
for j,v in enumerate(varyfactor):
    ax=plot_species(speciesNames[j],time,abundances[j],ax=axs,ls=Linestles[j],lw=3.0)


#ax2=ax.twinx()
#ax2.plot(time,temp)

#plot species returns the axis so we can further edit
axs.set(xscale='linear',yscale='log',xlim=(1,1e7),ylim=(5e-18,3e-3))
axs.set_ylabel('X/H',fontsize='x-large')
axs.set_xlabel('Time / Myears',fontsize='x-large')
axs.ticklabel_format(axis='x',useMathText=True)
axs.legend(loc='upper left',fontsize='xx-large',handlelength=20)
axs.set_title('Key for 5 levels of varied element',fontsize='xx-large')

#overwrite our previous plot
fig.savefig(plot_file)
