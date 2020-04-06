import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
# Specify the path for the ASCII data files
base_path = 'Benchmarking/'

# Specify the prefix and suffix for the model output data files
prefix = 'test'
suffix = '.out'


for radfield in ["low_rad/"]:
	fig,axes=plt.subplots(2,3,figsize=(16,9))
	axes=axes.flatten()
	for ax in axes:
		ax.set(xscale='log',yscale='log',xlim=(1e-4,1e2),xlabel="Av / mag",ylabel="$\Lambda$ / erg cm$^{-3}$ s^{-1}$")
	ax,ax2,ax3,ax4,ax5,ax6=axes

	path=base_path+radfield
	particle,*junk,av=np.loadtxt(path+"test.av.out",unpack=True,skiprows=1)
	temp=np.loadtxt(path+"test.temp.out",unpack=True,skiprows=1,usecols=[2])


	ax.plot(av,temp,label="UCLPDR")
	ax.set(ylabel="Gas Temperature / K")

	pdr_heat=pd.read_csv(f"{path}/test.abun.out",delimiter=" ")
	print(pdr_heat.head())
	ax4.plot(av,pdr_heat["X(H)"],label="H PDR",color="blue",ls="--")
	ax4.plot(av,pdr_heat["X(H2)"],label="H2 PDR",color="red",ls="--")
	ax4.set(ylabel="Abundance")

	#Cooling from UCL-PDR
	pdr_heat=pd.read_csv(f"{path}/test.cool.out",delimiter=" ")
	for col in pdr_heat.columns[1:]:
		if "named" not in col:
			if col=="Total":
				ax5.plot(av,pdr_heat[col],label=col,color="black")
			else:
				ax5.plot(av,pdr_heat[col],label=col)

	ax5.set(ylim=(1e-28,1e-18))
	ax5.legend()

	#heting from UCLPDR
	pdr_heat=pd.read_csv(f"{path}/test.heat.out",delimiter=" ")
	pdr_heat[pdr_heat.astype(float)<=0.0]=1.0e-50
	print(pdr_heat.head(5))
	for col in pdr_heat.columns[1:]:
		if "named" not in col:
			if col=="Total":
				ax6.plot(av,pdr_heat[col],label=col,color="black")
			else:
				ax6.plot(av,pdr_heat[col],label=col)

	ax6.set(ylim=(1e-28,1e-18))
	ax6.legend()

	

	model_no,cloud_size,av=np.loadtxt(path+"av_grid.dat",unpack=True,delimiter=",",skiprows=1)

	avs=[]
	hs=[]
	h2s=[]
	temps=[]
	cools=[]
	heats=[]
	for model in model_no:
		data=np.loadtxt(f"{path}av_grid/{model:.0f}.dat")
		avs.append(np.mean(data[-100:,3]))
		temp=np.mean(data[-100:,2])
		temps.append(temp)
		hs.append(np.mean(data[-100:,6]))
		h2s.append(np.mean(data[-100:,5]))


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



	ax.plot(avs,temps,label="UCLCHEM")
	ax.legend()
	
	cools=np.asarray(cools)
	for i in range(cools.shape[1]):
		ax2.plot(avs,cools[:,i],label=cool_labels[i])
	ax2.plot(avs,np.sum(cools,axis=1),label="Total",color="black")
	ax2.set(ylim=(1e-28,1e-18),xlim=(1e-4,1e2))
	ax2.legend()

	heats=np.asarray(heats)
	for i in range(heats.shape[1]):
		ax3.plot(avs,heats[:,i],label=heat_labels[i])
	ax3.plot(avs,np.sum(heats,axis=1),label="Total",color="black")		
	ax3.set(ylim=(1e-28,1e-18),xlim=(1e-4,1e2))
	ax3.legend()

	print(hs)
	ax4.plot(avs,hs,color="blue",label="H UCLCHEM")
	ax4.plot(avs,h2s,color="red",label="H2 UCLCHEM")
	ax4.legend()

	fig.savefig(path+"benchmark.png")