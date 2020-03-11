from __future__ import print_function
import uclchem
import numpy as np
import pandas as pd
import os
import time
import matplotlib.pyplot as plt


species_labels=["H2O","H2","H","CO","C+","O","CO+","CO2"]
models=[[100.0,10.0,1.0,2.0e6,"output/cold_diffuse.dat"],
		[1000000.0,10.0,1.0,2.0e6,"output/cold_dark.dat"],
		[100.0,100.0,1.0,2.0e6,"output/hot_diffuse.dat"],
		[1000000.0,100.0,1.0,2.0e6,"output/hot_dark.dat"]]

for model in models:
	in_dens,in_temp,uv,fin_time,file_out=model
	print(file_out)
	print(in_dens,in_temp,uv,fin_time,file_out)
	uclchem.temperaturetest(in_dens,in_temp,uv,fin_time,file_out)
	fig,ax=plt.subplots()
	time,dens,temp,*abundances=np.loadtxt(file_out,unpack=True)
	ax2=ax.twinx()
	for j,species in enumerate(species_labels):
		if species in ["H2","H","CO","C+","O"]:
			ax.plot(time,abundances[j],label=species)

	ax2.plot(time,temp,color="black",ls=":")
	ax.legend()
	ax.set(title="{0} K n$_H$={1} F$_{{UV}}$={2}".format(in_temp,in_dens,uv),yscale="log",ylim=(1e-16,1))
	fig.savefig(file_out.replace('.dat','.png'),dpi=300,bbox_inches='tight')