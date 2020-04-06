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

cloud_files={"square":"cloud-square.dat",
			"high_rad":"cloud-1d-n=1e3.dat",
			"low_rad":"cloud-1d-n=1e3.dat"}


for cloud_type in ["low_rad","high_rad"]:
	#load up density for each particle
	particle_df=pd.read_csv(f"{cloud_type}/test.prop.out",delimiter=" ")[["Particle","n_H"]]
	print(len(particle_df))
	print(particle_df.head())

	#then load h2 abundances and merge
	df=pd.read_csv(f"{cloud_type}/test.abun.out",delimiter=" ")[["Particle","X(H2)","X(C)"]]
	particle_df=particle_df.merge(df,on="Particle")
	print(len(particle_df))
	#then load ucl-pdr input file with distances
	df=pd.read_csv(cloud_files[cloud_type],delimiter=",")
	df=df[['Particle # ','      x (cm)      ']]
	df=df.rename({'Particle # ':"Particle",'      x (cm)      ':'x'},axis=1)
	particle_df=particle_df.merge(df,on="Particle")

	print(len(particle_df))
	#not sure why uclpdr uses the negative distance but we'll just use change in distance between particles
	particle_df["delta_x"]=0.0
	particle_df.loc[1:,"delta_x"]=particle_df.iloc[1:]["x"].values-particle_df.iloc[0:-1]["x"].values


	#can calulate change in column density between each particle as well
	particle_df["delta_h_col"]=particle_df["delta_x"].values*particle_df["X(H2)"].values*particle_df["n_H"].values
	particle_df["delta_c_col"]=particle_df["delta_x"].values*particle_df["X(C)"].values*particle_df["n_H"].values
	particle_df["delta_col_dens"]=particle_df["delta_x"].values*particle_df["n_H"].values


	#and then integrate
	particle_df["total_h_col"]=particle_df["delta_h_col"].cumsum()
	particle_df["total_c_col"]=particle_df["delta_c_col"].cumsum()
	particle_df["total_col_dens"]=particle_df["delta_col_dens"].cumsum()


	#integrate change in distance to get size
	particle_df["size"]=particle_df.delta_x.cumsum()


	#Then choose 48 evenly spaced points to sample
	step=int(np.floor(len(particle_df)/96))#only want 48 models
	print(step,len(particle_df))
	particle_df=particle_df.iloc[1:step*96+1:step] #only want 48 models
	particle_df.to_csv(f"{cloud_type}/model_tests.dat",index=False)

