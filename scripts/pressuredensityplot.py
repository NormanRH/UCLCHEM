#demonstration of plotfunctions. called from main UCLCHEM directory
#it reads full UCLCHEM output and saves a plot of the abudances of select species 

import numpy as np
from plotfunctions import *

#pick species, any number is fine
speciesNames=["H2S2"] #"H2","CO","H2O","CH3OH","#H2O","#CO","#CH3OH"]
input_file="../VaryFromSolar/outputfiles1/phase1-fullC-2.dat"
input_file2="../VaryFromSolar/outputfiles1/phase2-fullC-2.dat"
plot_file="../VaryFromSolar/KeyAndOtherPlots1/tempdensPlot.png"
input_file3="../VaryFromSolar/outputfiles1/static-fullS-2.dat"
#plot_file="../VaryFromSolar/KeyAndOtherPlots1/statictempdensPlot.png"


#call read_uclchem. 
t1,dens,temp,abundances=read_uclchem(input_file,speciesNames)
t2,dens2,temp2,abundances2=read_uclchem(input_file2,speciesNames)
t3,dens3,temp3,abundances3=read_uclchem(input_file3,speciesNames)

#write out to columnated output,
#write_cols("output/democolumns.dat",time,dens,abundances)

time = np.append(np.array(t1),np.array(t2) + t1[-1])
dens.extend(dens2)
temp.extend(temp2)

#plot species and save to test.png, alternatively send dens instead of time.
fig,axs=plt.subplots(figsize=(16,9)) 

#plot_species will make an axis if you don't provide one
#it will save if plotFile is given otherwise just returns figure and axis
#ax=plot_species(speciesNames,time,abundances,ax=axs)#,plotFile=plot_file)
axs.plot(time,dens,color='b',label='density',linewidth=3.0)
axs.plot(t3,dens3,color='tab:purple',label='static density',linewidth=2.0)

ax2=axs.twinx()
ax2.plot(time,temp,color='r',label='T',linewidth=3.0)
ax2.plot(t3,temp3,color='tab:orange',label='static T',linewidth=2.0)

#plot species returns the axis so we can further edit
axs.set(xscale='linear',yscale='log',xlim=(1,6.0e6),ylim=(1,1e6))
axs.set_xlabel('Time / Myears',fontsize='x-large')
axs.set_ylabel('n/cm^3',fontsize='x-large')
axs.set_title('Temperature/density Plot for collapse and static models',fontsize='xx-large')
axs.legend(loc='upper left',fontsize='large')
axs.ticklabel_format(axis='x',useMathText=True)#useMathText gives us a tidy 10^6 output
ax2.set(yscale='log',ylim=(1,600))
ax2.set_ylabel('T K',fontsize='x-large')
ax2.legend(loc='upper right',fontsize='large')
#overwrite our previous plot
fig.savefig(plot_file)
