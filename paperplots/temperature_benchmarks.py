import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.gridspec as gridspec

plt.style.use("thesis")
colors=sns.color_palette("deep",n_colors=6)


#####################################################################################################################

fig,axes=plt.subplots(3,2,figsize=(8.27,11.69))
axes=axes.flatten()

models={
	"low_rad":"1 Draine\n10$^3$ cm$^{-3}$",
	"high_cr":"1 Draine\n10$^3$ cm$^{-3}$\n100xCR",
	"10_1e3":"10 Draine\n10$^3$ cm$^{-3}$",
	"10_1e5.5":"10 Draine\n10$^{5.5}$ cm$^{-3}$",
	"1e5_1e3":"10$^5$ Draine\n10$^3$ cm$^{-3}$",
	"1e5_1e5.5":"10$^5$ Draine\n10$^{5.5}$ cm$^{-3}$"
}

for i,model in enumerate(models.keys()):
	uclchem=pd.read_csv(f"Benchmarking/uclchem-reduced/{model}.csv",usecols=["av","gasTemp","dustTemp"])
	uclpdr=pd.read_csv(f"Benchmarking/uclpdr/{model}_output.csv",usecols=["Av","T_g","T_d"])

	axes[i].plot(uclchem["av"],uclchem["gasTemp"],label="T$_g$ UCLCHEM",color="black")
	axes[i].plot(uclpdr["Av"],uclpdr["T_g"],label="T$_g$ UCL-PDR",ls="--",color="black")
	axes[i].plot(uclchem["av"],uclchem["dustTemp"],label="T$_d$ UCLCHEM",color="orange")
	axes[i].plot(uclpdr["Av"],uclpdr["T_d"],label="T$_d$ UCL-PDR",ls="--",color="orange")

	axes[i].text(0.99,0.99,models[model],horizontalalignment="right",verticalalignment="top",transform=axes[i].transAxes)

	axes[i].set(xscale='log',xlabel="A$_v$ / Mag", ylabel="Temperature / K")


axes[0].legend()
fig.subplots_adjust(wspace=0.3)
fig.savefig("paperplots/temperature-grid.pdf",type="PDF",bbox_inches="tight",dpi=300)

fig,ax=plt.subplots()
model="sine"


#####################################################################################################################

fig = plt.figure(tight_layout=True,figsize=(4.1,5.5))
gs = gridspec.GridSpec(3, 1)

ax = fig.add_subplot(gs[0, 0])

uclchem=pd.read_csv(f"Benchmarking/uclchem-reduced/{model}.csv")
uclpdr=pd.read_csv(f"Benchmarking/uclpdr/{model}_output.csv")
ax.plot(uclchem["av"],uclchem["gasTemp"],label="T$_g$ UCLCHEM",color="black")
ax.plot(uclpdr["Av"],uclpdr["T_g"],label="T$_g$ UCL-PDR",ls="--",color="black")
ax.plot(uclchem["av"],uclchem["dustTemp"],label="T$_d$ UCLCHEM",color="orange")
ax.plot(uclpdr["Av"],uclpdr["T_d"],label="T$_d$ UCL-PDR",ls="--",color="orange")
ax.set(xscale='log',xlabel="A$_v$ / Mag", ylabel="Temperature / K")

axis = fig.add_subplot(gs[1:,0])


uclpdr_species=["H_abun","H2_abun","H+_abun","C+_abun","C_abun","e-_abun"]
for i,col in enumerate(uclpdr_species):
	axis.plot(uclpdr["Av"],uclpdr[col],color=colors[i],ls="--",label="")

uclchem_species=["H","H2","H+","C+","C","E-"]
for i,col in enumerate(uclchem_species):
	axis.plot(uclchem["av"],uclchem[col],color=colors[i],label=col)
axis.set(xscale='log',yscale='log',xlabel="A$_v$ / mag",ylabel="Abundance")
axis.legend(loc=2)

fig.savefig("paperplots/sinewave.pdf",type="PDF",bbox_inches="tight")