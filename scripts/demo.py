#demonstration of plotfunctions. called from main UCLCHEM directory
#it reads full UCLCHEM output and saves a plot of the abudances of select species 

from plotfunctions import *
import pandas as pd

#pick species, any number is fine
#zspeciesNames=["CO","CS","H2O","#CO","#CH3OH","NH3"]
speciesNames=["H2O","H2","H","CO","C+","O","CO+","CO2"]
file_name="column"
plot_file="output/"+file_name+".png"
#call read_uclchem. 
time,dens,temp,*abundances=np.loadtxt("output/"+file_name+".dat",unpack=True)
temp=pd.Series(temp).rolling(window=15).mean().values
#write out to columnated output,
#write_cols("output/democolumns.dat",time,dens,abundances)

#plot species and save to test.png, alternatively send dens instead of time.
axis,fig=plot_species(speciesNames,time,abundances,plotFile=plot_file)
ax2=axis.twinx()
ax2.plot(time,temp,ls=":",color="black")
#plot species returns the axis so we can further edit
axis.set(ylim=(1e-20,1e-3),xscale='log')
axis.set_title('Coolants')
fig.savefig(plot_file)
