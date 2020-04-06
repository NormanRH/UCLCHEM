import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
# Specify the path for the ASCII data files
base_path = 'Benchmarking/'

# Specify the prefix and suffix for the model output data files
prefix = 'test'
suffix = '.out'


for radfield in ["low_rad/"]:#,"high_rad/","low_rad/"]:
	fig,axes=plt.subplots(2,3,figsize=(16,9))
	axes=axes.flatten()
	for ax in axes:
		ax.set(xscale='log',yscale='log',xlabel="Av / mag",ylabel="$\Lambda$ / erg cm$^{-3}$ $s^{-1}$")
	ax,ax2,ax3,ax4,ax5,ax6=axes

	path=base_path+radfield
	particle,*junk,av=np.loadtxt(path+"test.av.out",unpack=True,skiprows=1)
	temp,dusttemp=np.loadtxt(path+"test.prop.out",unpack=True,skiprows=1,usecols=[3,4])


	ax.plot(av,temp,label="T$_{g}$ PDR",color="blue",ls="--")
	ax.plot(av,dusttemp,label="T$_{d}$ PDR",color="orange",ls="--")
	ax.set(ylabel="Gas Temperature / K")

	pdr_heat=pd.read_csv(f"{path}/test.abun.out",delimiter=" ")
	print(pdr_heat.head())
	ax4.plot(av,pdr_heat["X(H)"],label="H PDR",ls="--")
	ax4.plot(av,pdr_heat["X(H2)"],label="H2 PDR",ls="--")
	ax4.plot(av,pdr_heat["X(H+)"],label="H+ PDR",ls="--")
	#ax4.plot(av,pdr_heat["X(C+)"],label="C+ PDR,ls="--")
	#ax4.plot(av,pdr_heat["X(O)"],label="O PDR",ls="--")
	ax4.plot(av,pdr_heat["X(e-)"],label="E- PDR",ls="--")
	ax4.set(ylabel="Abundance")

	#Cooling from UCL-PDR
	pdr_heat=pd.read_csv(f"{path}/test.cool.out",delimiter=" ")
	for col in pdr_heat.columns[1:]:
		if "named" not in col:
			if col=="Total":
				ax5.plot(av,pdr_heat[col],label=col,color="black")
			else:
				ax5.plot(av,pdr_heat[col],label=col)

	ax5.text(0.05,0.95,"UCLPDR",transform=ax5.transAxes)
	#ax5.legend()

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

	ax6.text(0.05,0.95,"UCLPDR",transform=ax6.transAxes)
	#ax6.legend()

	

	model_no,cloud_size,av=np.loadtxt(path+"av_grid.dat",unpack=True,delimiter=",",skiprows=1,usecols=[0,1,2])
	species=["H","H2","E-","H+"]
	abundances=[[] for  spec in species]
	avs=[]
	temps=[]
	dusttemps=[]
	cools=[]
	heats=[]
	for model in model_no:
		data=pd.read_csv(f"{path}av_grid/{model:.0f}.dat").rename(str.strip,axis=1)
		avs.append(av[int(model)])
		temp=np.mean(data["gasTemp"].values[-100:])
		temps.append(temp)
		#dust
		temp=np.mean(data["dustTemp"].values[-100:])
		dusttemps.append(temp)

		#"H2O,H2,H,H+,HE+,C,C+,O,O+,CO,CO+,CO2,#CO,#H2O,E-"
		for i,spec in enumerate(species):
			abundances[i].append(np.mean(data[spec].values[-100:]))



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


	print(temps,avs)


	ax.plot(avs,temps,label="T$_{g}$ UCLCHEM",color="blue")
	ax.plot(avs,dusttemps,label="T$_{d}$ UCLCHEM",color="orange")
	ax.legend()
	
	cools=np.asarray(cools)
	for i in range(cools.shape[1]):
		ax2.plot(avs,cools[:,i],label=cool_labels[i])
	ax2.plot(avs,np.sum(cools,axis=1),label="Total",color="black")
	ax2.text(0.05,0.95,"UCLCHEM",transform=ax2.transAxes)
	ax2.legend()


	heats=np.asarray(heats)
	for i in range(heats.shape[1]):
		ax3.plot(avs,heats[:,i],label=heat_labels[i])
	ax3.plot(avs,np.sum(heats,axis=1),label="Total",color="black")		
	ax3.set()
	ax3.text(0.05,0.95,"UCLCHEM",transform=ax3.transAxes)

	ax3.legend()

	for i,spec in enumerate(species):
		label=spec+" UCLCHEM"
		ax4.plot(avs,abundances[i],label=label)
	ax4.legend()

	for axis in [ax2,ax3,ax5,ax6]:
		axis.set(xlim=(np.min(avs),np.max(avs)),ylim=(1e-28,1e-21))

	fig.savefig(path+"benchmark.png",dpi=300,bbox_inches="tight")