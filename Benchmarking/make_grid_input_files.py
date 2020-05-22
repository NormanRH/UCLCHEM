##################################################################################
# My wrap.f90 takes h2col and coldens as inputs
# This script creates a file with those values to be fed into h2col_grid.py
# allowing an equivalent model to a given UCLPDR model to be run.
#
# requires the input particle file for the UCLPDR model and some out of the outputs
##################################################################################

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

cloud_files={"square":"square.cloud",
			"sine":"sine.cloud",
			"10_linear_decrease":"linear-density-decrease.cloud",
			"10_linear_increase":"linear-density-increase.cloud",
			"high_rad":"1e3.cloud",
			"low_rad":"1e3.cloud",
			"low_rad_fixed":"1e3.cloud",
			"high_cr":"1e3.cloud",
			"10_1e3":"1e3.cloud",
			"10_1e5":"1e5.5.cloud",
			"1e5_1e3":"1e3.cloud",
			"1e5_1e5":"1e5.5.cloud",
			"low_metallicity":"1e3.cloud",
			"fixed_cooling":"1e5.5.cloud"}



for cloud_type in ["low_metallicity"]:#square","10_linear_increase","10_linear_decrease","10_1e5","1e5_1e5","10_1e3","1e5_1e3","high_cr
	#load up density for each particle
	particle_df=pd.read_csv(f"Benchmarking/uclpdr/{cloud_type}_output.csv")[["Particle","T_g","n_H","H2_abun","C_abun","FUV"]]

	#get surface uv for all particles
	particle_df["FUV"]=particle_df.loc[1,"FUV"]

	#then load ucl-pdr input file with distances
	df=pd.read_csv(f"Benchmarking/cloud_files/{cloud_files[cloud_type]}",delimiter=",")[["Particle#","x(cm)"]]
	df=df.rename({"Particle#":"Particle","x(cm)":"x"},axis=1)
	particle_df=particle_df.merge(df,on="Particle")

	#not sure why uclpdr uses the negative distance but we'll just use change in distance between particles
	particle_df["delta_x"]=0.0
	particle_df.loc[1:,"delta_x"]=particle_df.iloc[1:]["x"].values-particle_df.iloc[0:-1]["x"].values
	print(particle_df.head())

	#can calulate change in column density between each particle as well
	particle_df["delta_h2_col"]=particle_df["delta_x"].values*particle_df["H2_abun"].values*particle_df["n_H"].values
	particle_df["delta_c_col"]=particle_df["delta_x"].values*particle_df["C_abun"].values*particle_df["n_H"].values
	particle_df["delta_col_dens"]=particle_df["delta_x"].values*particle_df["n_H"].values


	#and then integrate
	particle_df["total_h2_col"]=particle_df["delta_h2_col"].cumsum()
	particle_df["total_c_col"]=particle_df["delta_c_col"].cumsum()
	particle_df["total_col_dens"]=particle_df["delta_col_dens"].cumsum()


	#integrate change in distance to get size
	particle_df["size"]=particle_df.delta_x.cumsum()

	n_samples=96
	#Then choose 48 evenly spaced points to sample
	step=int(np.floor(len(particle_df)/n_samples))#only want 48 models
	print(step,len(particle_df))
	particle_df=particle_df.iloc[0:step*n_samples+1:step] #only want 48 models
	particle_df.to_csv(f"Benchmarking/grid_inputs/{cloud_type}.csv",index=False)

