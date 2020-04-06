import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
# Specify the path for the ASCII data files
base_path = 'Benchmarking/'

# Specify the prefix and suffix for the model output data files
prefix = 'test'
suffix = '.out'


for radfield in ["low_rad/"]:#,"high_rad/","low_rad/"]:
	path=base_path+radfield

	fig,axes=plt.subplots(1,2,figsize=(16,9))
	axes=axes.flatten()
	for ax in axes:
		ax.set(xscale='log',yscale='log',xlim=(1e-4,1e2),xlabel="Av / mag",ylabel="$\Lambda$ / erg cm$^{-3}$ $s^{-1}$")
	ax,ax2=axes
	particle,*junk,av=np.loadtxt(path+"test.av.out",unpack=True,skiprows=1)

	#Cooling from UCL-PDR
	pdr_heat=pd.read_csv(f"{path}/test.cool.out",delimiter=" ")
	for col in pdr_heat.columns[1:]:
		if "named" not in col:
			if col=="Total":
				ax.plot(av,pdr_heat[col],color="black",ls="--")
			else:
				ax.plot(av,pdr_heat[col],ls="--")

	ax.set(ylim=(1e-26,1e-20))
	#ax5.legend()

	#heting from UCLPDR
	pdr_heat=pd.read_csv(f"{path}/test.heat.out",delimiter=" ")
	pdr_heat[pdr_heat.astype(float)<=0.0]=1.0e-50
	print(pdr_heat.head(5))
	for col in pdr_heat.columns[1:]:
		if "named" not in col:
			if col=="Total":
				ax2.plot(av,pdr_heat[col].values,color="black",ls="--")
			else:
				ax2.plot(av,pdr_heat[col].values,ls="--")

	ax2.set(ylim=(1e-26,1e-20))
	#ax6.legend()

	

	model_no,cloud_size,av=np.loadtxt(path+"av_grid.dat",unpack=True,delimiter=",",skiprows=1,usecols=[0,1,2])

	avs=[]
	hs=[]
	h2s=[]
	temps=[]
	dusttemps=[]
	cools=[]
	heats=[]
	for model in model_no:
		data=np.loadtxt(f"{path}av_grid/{model:.0f}.dat")
		avs.append(av[int(model)])
		#cooling
		cool_labels=["Ly-alpha","C+","O","C","CO","p-H2","o-H2"]
		data=np.loadtxt(f"{path}av_grid/cooling/{model:.0f}.dat",skiprows=1)
		cool=np.median(data[-100:,:],axis=0)
		cool[np.where(cool==0.0)]==1.0e-50
		cools.append(np.median(data[-100:,:],axis=0))
		
		#heating
		heat_labels=["Photoelectric","H2Formation","FUVPumping","Photodissociation","Cionization","CRheating","Chemheating","turbHeating","gasGrainColls"]
		data=np.loadtxt(f"{path}/av_grid/heating/{model:.0f}.dat",skiprows=1)
		heats.append(np.median(data[-100:,:],axis=0))

	cools=np.asarray(cools)
	for i in range(cools.shape[1]):
		ax.plot(avs,cools[:,i],label=cool_labels[i])
	ax.plot(avs,np.sum(cools,axis=1),label="Total",color="black")
	ax.set(ylim=(1e-28,1e-21),xlim=(1e-4,1e2))
	ax.legend()


	heats=np.asarray(heats)
	for i in range(heats.shape[1]):
		ax2.plot(avs,heats[:,i],label=heat_labels[i])
	ax2.plot(avs,np.sum(heats,axis=1),label="Total",color="black")		
	ax2.set(ylim=(1e-28,1e-20),xlim=(1e-4,1e2))
	ax2.legend()
	fig.savefig(path+"heating-cooling.png",dpi=300,bbox_inches="tight")