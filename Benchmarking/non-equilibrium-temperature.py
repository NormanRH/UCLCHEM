from glob import glob
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np


files=glob("Benchmarking/*/av_grid/*.dat")

data=pd.read_csv(files[0]).rename(str.strip,axis=1)[["Time","gasTemp"]]
data["gasTemp"]=data["gasTemp"].values/np.median(data["gasTemp"].values[-100:])
	
for i,file in enumerate(files[1:]):
	temp=pd.read_csv(file).rename(str.strip,axis=1)[["Time","gasTemp"]]
	temp["gasTemp"]=temp["gasTemp"].values/np.median(temp["gasTemp"].values[-100:])
	data=data.append(temp)


fig,ax=plt.subplots(figsize=(16,9))
sns.lineplot(data=data,x="Time",y="gasTemp")
ax.set(xscale='log',xlabel="Time / year",ylabel="Temperature/Final Temperature",xlim=(1,2e6))
fig.savefig("Benchmarking/equilibrium-temps.png")