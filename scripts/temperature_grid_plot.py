import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

models=pd.read_csv("output/grid.dat",delimiter=" ")
models["dens"]=models["dens"].map(np.log10)
max_temp=[]
avg_temp=[]
start_temp=[]
start_dens=[]
final_temp=[]
passes=[]

for i,row in models.iterrows():
    try:
        data=np.loadtxt(row["file"])
        avg_temp.append(data[:,2].mean())
        final_temp.append(data[-1,2])
        data[:,2]=data[:,2]-row["temp"]
        passes.append(True)
        max_temp.append(data[np.argmax(np.abs(data[:,2])),2].max()) #time,dens,temp, abunds
    except:
        print(row["file"]," failed")
        print(row)
        passes.append(False)
models=models[passes]

models["Max Temp Diff"]=max_temp
models["Average Temp"]=avg_temp
models["Final Temp"]=final_temp


plots=["Final Temp","Average Temp"]
fig,ax=plt.subplots(len(plots),3,figsize=(16,9))

cmap=sns.cubehelix_palette(start=1,rot=0.03,gamma=1.36,hue=2,as_cmap=True,reverse=True)
for i,plot in enumerate(plots):
    for j, uv in enumerate(models["uv"].unique()):
        average_map=pd.pivot_table(models[models["uv"]==uv],index="temp",columns="dens",values=plot).fillna(0.0)
        average_map=average_map.sort_index(ascending=True)

        #heatmap = ax[i,j].pcolormesh(average_map, alpha=1)
        heatmap=sns.heatmap(average_map,ax=ax[i,j],cbar_kws={'label': plot},cmap=cmap)
        ax[i,j].invert_yaxis()
        ax[i,j].set_yticks(np.arange(average_map.shape[0]) + 0.5, minor=False)
        ax[i,j].set_xticks(np.arange(average_map.shape[1]) + 0.5, minor=False)
        ax[i,j].set_xticklabels(average_map.columns, minor=False)
        ax[i,j].set_yticklabels(average_map.index, minor=False)

        ax[i,j].set(ylabel="Initial Temperature / K",xlabel="log(n$_H$ / cm$^{-3}$)",title="UV = {0} Habing".format(uv))
plt.subplots_adjust(wspace=0.3,hspace=0.3)
plt.show()
fig.savefig("output/temperature_grid.png",dpi=300,bbox_inches='tight')