import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
# Specify the path for the ASCII data files
uclchem_species=["H","H2","H+","C+","E-"]
uclchem_cools=["Ly-alpha","C+ ","O ","C ","CO ","p-H2","o-H2"]
uclchem_heats=["Photoelectric","H2Formation","FUVPumping","Photodissociation","Cionization","CRheating","Chemheating","turbHeating","gasGrainColls"]
colors=sns.color_palette("deep",n_colors=len(uclchem_heats)+1)

def pdr_heating_plot(axis,av_points,heat_file,linestyle="-"):
	pdr_heat=pd.read_csv(heat_file,delimiter=" ")
	pdr_heat[pdr_heat.astype(float)<=0.0]=1.0e-50
	i=0
	for col in pdr_heat.columns[1:]:
		if "named" not in col:
			if col=="Total":
				axis.plot(av_points,pdr_heat[col],label=col,color="black",ls=linestyle)
			else:
				axis.plot(av_points,pdr_heat[col],label=col,color=colors[i],ls=linestyle)
				i=i+1
	axis.set(xscale='log',yscale='log',xlabel="Av / mag",ylabel="$\Lambda$ / erg cm$^{-3}$ $s^{-1}$")
	axis.set_ylim(bottom=1e-32)
	return axis


def pdr_abundance_plot(axis,av_points,abund_file):
	abund_df=pd.read_csv(abund_file,delimiter=" ")

	axis.plot(av_points,abund_df["X(H)"],color=colors[0],label="H PDR",ls="--")
	axis.plot(av_points,abund_df["X(H2)"],color=colors[1],label="H2 PDR",ls="--")
	axis.plot(av_points,abund_df["X(H+)"],color=colors[2],label="H+ PDR",ls="--")
	axis.plot(av_points,abund_df["X(C+)"],color=colors[3],label="C+ PDR",ls="--")
		#ax4.plot(av,pdr_heat["X(O)"],label="O PDR",ls="--")
	axis.plot(av_points,abund_df["X(e-)"],color=colors[4],label="E- PDR",ls="--")
	axis.set(xscale='log',yscale='log',xlabel="Av / mag",ylabel="Abundance")
	return axis

def create_uclchem_df(path):
	model_no,cloud_size,av=np.loadtxt(path+"/av_grid.dat",unpack=True,delimiter=",",skiprows=1,usecols=[0,1,2])
	species=["H","H2","H+","C+","E-"]
	abundances=[[] for  spec in species]
	avs=[]
	temps=[]
	dusttemps=[]
	cools=[]
	heats=[]
	for model in model_no:
		data=pd.read_csv(f"{path}/av_grid/{model:.0f}.dat").rename(str.strip,axis=1)
		avs.append(av[int(model)])
		temp=np.median(data["gasTemp"].values[-100:])
		temps.append(temp)
		#dust
		temp=np.median(data["dustTemp"].values[-100:])
		dusttemps.append(temp)

		#"H2O,H2,H,H+,HE+,C,C+,O,O+,CO,CO+,CO2,#CO,#H2O,E-"
		for i,spec in enumerate(species):
			abundances[i].append(np.mean(data[spec].values[-100:]))


		#cooling
		if "fixed_cooling" in path:
			cool=[0.0,0.0,0.0,0.0,0.0,0.0,0.0]
		else:
			data=np.loadtxt(f"{path}/av_grid/cooling/{model:.0f}.dat",skiprows=1)
			cool=np.median(data[-100:,:],axis=0)
			cool[np.where(cool<1.0e-50)]=1.0e-50
		cools.append(cool)
		
		#heating
		data=np.loadtxt(f"{path}/av_grid/heating/{model:.0f}.dat",skiprows=1)
		heats.append(np.median(data[-100:,:],axis=0))

	cools=np.asarray(cools)
	heats=np.asarray(heats)
	#create a dictionary of each variable and list of particle values
	uclchem_df={"av":avs,"gasTemp":temps,"dustTemp":dusttemps}
	for i in range(cools.shape[1]):
		uclchem_df[uclchem_cools[i]]=cools[:,i]
	for i in range(heats.shape[1]):
		uclchem_df[uclchem_heats[i]]=heats[:,i]
	for i,spec in enumerate(species):
		uclchem_df[spec]=abundances[i]

	uclchem_df=pd.DataFrame(uclchem_df)
	return uclchem_df

def uclchem_heating_plot(uclchem_df,axis,type):
	if type=="cooling":
		columns=uclchem_cools
	elif type=="heating":
		columns=uclchem_heats
	else:
		print("unknown plot type for uclchem!")
		exit
	for i,col in enumerate(columns):
		axis.plot(uclchem_df["av"],uclchem_df[col],color=colors[i])
	axis.plot(uclchem_df["av"],uclchem_df[columns].sum(axis=1),label="Total",color="black")
	axis.set(xscale='log',yscale='log',xlabel="Av / mag",ylabel="$\Lambda$ / erg cm$^{-3}$ $s^{-1}$")
	axis.set_ylim(bottom=1e-32)

	return axis

def uclchem_abund_plot(uclchem_df,axis):
	columns=uclchem_species
	for i,col in enumerate(columns):
		axis.plot(uclchem_df["av"],uclchem_df[col],color=colors[i])
	axis.set(xscale='log',yscale='log',xlabel="Av / mag",ylabel="Abundance")
	return axis

def temperature_plot(axis,pdr_av,pdr_temp_file,uclchem_df,legend=True):
	temp,dusttemp=np.loadtxt(pdr_temp_file,unpack=True,skiprows=1,usecols=[3,4],delimiter=",")
	axis.plot(pdr_av,temp,ls="--",label="PDR T$_{g}$",color=colors[0])
	axis.plot(pdr_av,dusttemp,ls="--",label="PDR T$_{d}$",color=colors[3])
	axis.plot(uclchem_df["av"],uclchem_df["gasTemp"],label="UCLCHEM T$_{g}$",color=colors[0])
	axis.plot(uclchem_df["av"],uclchem_df["dustTemp"],label="UCLCHEM T$_{d}$",color=colors[3])
	axis.set(xscale='log',ylabel="Temperature / K",xlabel="Av / mag")
	if legend:
		axis.legend()
	return axis
