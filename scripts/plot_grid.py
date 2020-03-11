import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

models=pd.read_csv("output/grid.dat",delimiter=" ")
outputFolder="output/grid_plots/"


species_labels=["H2O","H2","H","CO","C+","O","CO+","CO2"]#,"#CO","#H2O"]
for i,row in models.iterrows():
    print(row["file"])
    try:
        time,dens,temp,*abundances=np.loadtxt(row["file"],unpack=True)
        read=True
    except:
        read=False
        print(row["file"],"failed")
    if read:
        fig,ax=plt.subplots()
        ax2=ax.twinx()
        for j,species in enumerate(species_labels):
            if species in ["H2","H","CO","C+","O"]:
                ax.plot(time,abundances[j],label=species)
        
        ax2.plot(time,temp,color="black",ls=":")
        ax.legend()
        ax.set(title="{0} K n$_H$={1} F$_{{UV}}$={2}".format(row["temp"],row["dens"],row["uv"]),yscale="log",ylim=(1e-16,1))
        fig.savefig(outputFolder+str(row["model"]),dpi=300,bbox_inches='tight')