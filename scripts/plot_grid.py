import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

models=pd.read_csv("output/av_grid.dat",delimiter=",")
outputFolder="output/grid_plots/"


species_labels=["H2O","H2","H","CO","C+","O","CO+","CO2","#CO","#H2O"]
for i,row in models.iterrows():
    print(row["model"])
    time,dens,temp,av,*abundances=np.loadtxt(f"output/av_grid/{row.model:.0f}.dat",unpack=True)
    fig,ax=plt.subplots(figsize=(16,9))
    ax2=ax.twinx()
    for j,species in enumerate(species_labels):
        if species in ["H2","H","CO","C+","O"]:
            ax.plot(time,abundances[j],label=species)
    
    ax2.plot(time,temp,color="black",ls=":")
    ax.legend()
    ax.set(title=f"{row.av} Mags",yscale="log",ylim=(1e-16,1))
    fig.savefig(outputFolder+str(row["model"])+".png",dpi=300,bbox_inches='tight')