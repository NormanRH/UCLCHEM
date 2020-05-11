from plotting_functions import *

base_path = 'Benchmarking/'

for radfield in ["square"]:#"10_1e5.5","1e5_1e3","1e5_1e5.5","high_cr","low_rad_fixed","10_1e3",
	#set file variables
	path=base_path+radfield
	pdr_heat_file=f"{path}/{radfield}.heat.out"
	pdr_cool_file=f"{path}/{radfield}.cool.out"
	pdr_abund_file=f"{path}/{radfield}.abun.out"
	pdr_temp_file=f"{path}/{radfield}.prop.out"

	#load av at each particle so we can plot av on x axis.
	particle,*junk,pdr_av_points=np.loadtxt(f"{path}/{radfield}.av.out",unpack=True,skiprows=1,delimiter=",")
	#uclchem_df=create_uclchem_df(path)
	#uclchem_df.to_csv(f"{path}/uclchem_outputs.dat")
	uclchem_df=pd.read_csv(f"{path}/uclchem_outputs.dat")

	#single cooling plot
	fig,axis=plt.subplots(figsize=(16,9))
	axis=pdr_heating_plot(axis,pdr_av_points,pdr_cool_file,linestyle="--")
	axis.legend(fontsize="small")
	axis=uclchem_heating_plot(uclchem_df,axis,"cooling")
	fig.savefig(path+f"/{radfield}-cooling.png",bbox_inches="tight",dpi=300)

	#single heating plot
	fig,axis=plt.subplots(figsize=(16,9))
	axis=pdr_heating_plot(axis,pdr_av_points,pdr_heat_file,linestyle="--")
	axis.legend(fontsize="small")
	axis=uclchem_heating_plot(uclchem_df,axis,"heating")
	fig.savefig(path+f"/{radfield}-heating.png",bbox_inches="tight",dpi=300)

	#single temperature plot
	fig,axis=plt.subplots(figsize=(16,9))
	axis=temperature_plot(axis,pdr_av_points,pdr_temp_file,uclchem_df)
	fig.savefig(path+f"/{radfield}-temperature.png",bbox_inches="tight",dpi=300)




	fig,axes=plt.subplots(2,3,figsize=(16,9))
	axes=axes.flatten()
	
	axes[0]=temperature_plot(axes[0],pdr_av_points,pdr_temp_file,uclchem_df)
	axes[1]=uclchem_heating_plot(uclchem_df,axes[1],"cooling")
	axes[1].legend(fontsize="small")
	axes[2]=uclchem_heating_plot(uclchem_df,axes[2],"heating")
	axes[2].legend(fontsize="small")
	axes[3]=uclchem_abund_plot(uclchem_df,axes[3])
	axes[3].legend(fontsize="small")
	axes[3]=pdr_abundance_plot(axes[3],pdr_av_points,pdr_abund_file)
	axes[4]=pdr_heating_plot(axes[4],pdr_av_points,pdr_cool_file)
	axes[5]=pdr_heating_plot(axes[5],pdr_av_points,pdr_heat_file)

	axes[1].text(0.01,0.95,"UCLCHEM",transform=axes[1].transAxes)
	axes[2].text(0.01,0.95,"UCLCHEM",transform=axes[2].transAxes)
	axes[4].text(0.01,0.95,"UCL-PDR",transform=axes[4].transAxes)
	axes[5].text(0.01,0.95,"UCL-PDR",transform=axes[5].transAxes)


	fig.savefig(path+f"/{radfield}-benchmark.png",dpi=300,bbox_inches="tight")