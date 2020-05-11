from plotting_functions import *

base_path = 'Benchmarking/'

for radfield in ["10_1e5.5","1e5_1e3","1e5_1e5.5","high_cr","10_1e3","low_rad"]:
	
	uclpdr_df=pd.read_csv(f"{base_path}uclpdr/{radfield}_output.csv")
	uclchem_df=pd.read_csv(f"{base_path}uclchem/{radfield}.csv")

	# #single cooling plot
	# fig,axis=plt.subplots(figsize=(16,9))
	# axis=pdr_heating_plot(axis,uclpdr_df,subscript="cool",linestyle="--")
	# axis.legend(fontsize="small")
	# axis=uclchem_heating_plot(uclchem_df,axis,"cooling")
	# fig.savefig(path+f"/{radfield}-cooling.png",bbox_inches="tight",dpi=300)

	# #single heating plot
	# fig,axis=plt.subplots(figsize=(16,9))
	# axis=pdr_heating_plot(axis,uclpdr_df,subscript="heat",linestyle="--")
	# axis.legend(fontsize="small")
	# axis=uclchem_heating_plot(uclchem_df,axis,"heating")
	# fig.savefig(path+f"/{radfield}-heating.png",bbox_inches="tight",dpi=300)

	# #single temperature plot
	# fig,axis=plt.subplots(figsize=(16,9))
	# axis=temperature_plot(axis,uclpdr_df,uclchem_df)
	# fig.savefig(path+f"/{radfield}-temperature.png",bbox_inches="tight",dpi=300)



	fig,axes=plt.subplots(2,3,figsize=(16,9))
	axes=axes.flatten()
	
	axes[0]=temperature_plot(axes[0],uclpdr_df,uclchem_df)
	axes[1]=uclchem_heating_plot(uclchem_df,axes[1],"cooling")
	axes[1].legend(fontsize="small")
	axes[2]=uclchem_heating_plot(uclchem_df,axes[2],"heating")
	axes[2].legend(fontsize="small")
	axes[3]=uclchem_abund_plot(uclchem_df,axes[3])
	axes[3].legend(fontsize="small")
	axes[3]=pdr_abundance_plot(axes[3],uclpdr_df)
	axes[4]=pdr_heating_plot(axes[4],uclpdr_df,subscript="cool",linestyle="--")
	axes[5]=pdr_heating_plot(axes[5],uclpdr_df,subscript="heat",linestyle="--")

	axes[1].text(0.01,0.95,"UCLCHEM",transform=axes[1].transAxes)
	axes[2].text(0.01,0.95,"UCLCHEM",transform=axes[2].transAxes)
	axes[4].text(0.01,0.95,"UCL-PDR",transform=axes[4].transAxes)
	axes[5].text(0.01,0.95,"UCL-PDR",transform=axes[5].transAxes)


	fig.savefig(f"{base_path}comparison_plots/{radfield}.png",dpi=300,bbox_inches="tight")