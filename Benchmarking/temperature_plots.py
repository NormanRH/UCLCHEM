from plotting_functions import *

base_path = 'Benchmarking/'

fig,axes=plt.subplots(3,2,figsize=(8.27,11.7))
axes=axes.flatten()

labels={"10_1e3":"10 Draine \n 10$^3$ cm$^{-3}$",
"10_1e5.5":"10 Draine \n 10$^5.5$ cm$^{-3}$",
"1e5_1e3":"10$^5$ Draine \n 10$^3$ cm$^{-3}$",
"1e5_1e5.5":"10$^5$ Draine \n 10$^5.5$ cm$^{-3}$",
"high_cr":"1 Draine \n 10$^3$ cm$^{-3}$\n 100xCR",
"low_rad":"1 Draine \n 10$^3$ cm$^{-3}$"}


for i,radfield in enumerate(["low_rad","high_cr","10_1e3","10_1e5.5","1e5_1e3","1e5_1e5.5",]):
	axis=axes[i]
	#set file variables
	path=base_path+radfield
	pdr_heat_file=f"{path}/{radfield}.heat.out"
	pdr_cool_file=f"{path}/{radfield}.cool.out"
	pdr_abund_file=f"{path}/{radfield}.abun.out"
	pdr_temp_file=f"{path}/{radfield}.prop.out"

	#load av at each particle so we can plot av on x axis.
	particle,*junk,pdr_av_points=np.loadtxt(f"{path}/{radfield}.av.out",unpack=True,skiprows=1)
	#uclchem_df=create_uclchem_df(path)
	#uclchem_df.to_csv(f"{path}/uclchem_outputs.dat")
	uclchem_df=pd.read_csv(f"{path}/uclchem_outputs.dat")


	axis=temperature_plot(axis,pdr_av_points,pdr_temp_file,uclchem_df,legend=False)
	if i==0:
		axis.legend(fontsize="small")

	axis.text(0.95,0.95,labels[radfield],transform=axis.transAxes,horizontalalignment="right",verticalalignment="top")



fig.savefig("temperature_plots.pdf",type="pdf",dpi=300,bbox_inches="tight")