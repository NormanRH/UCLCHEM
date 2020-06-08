import uclchem
import numpy as np
import time
from multiprocessing import Pool
import pandas as pd
from os import path,mkdir
import random

def run_model(a):
	i,params=a
	outSpecies="H2O,H2,H,H+,HE+,C,C+,O,O+,CO,CO+,E-"
	uclchem.general(params,outSpecies,f"Dimensionality/raw/",f"{i:.0f}")
	return 1

def generate_elements(param_dict):
	min_scale=0.001
	max_scale=1.0
	param_dict["fc"]=random.uniform(min_scale,max_scale)*2.6e-4
	param_dict["fo"]=random.uniform(min_scale,max_scale)*4.6e-4
	param_dict["fn"]=random.uniform(min_scale,max_scale)*6.1e-5
	param_dict["fs"]=random.uniform(min_scale,max_scale)*1.3e-5
	param_dict["fsi"]=random.uniform(min_scale,max_scale)*1.0e-7
	param_dict["fcl"]=random.uniform(min_scale,max_scale)*3.162e-7
	param_dict["fmg"]=random.uniform(min_scale,max_scale)*3.98e-5
	return param_dict

def generate_column_dens(param_dict):
	col_dens=random.uniform(0.0,24.0)
	col_dens=10.0**col_dens

	#after a certain av, c col is high.
	if col_dens>1e22:
		c_col=(10.0**random.uniform(16.7,19.7))
		ion=0
	else:
		c_col=(10.0**random.uniform(-15.0,-4.0))
		if c_col>1e-6:
			ion=1
		else:
			ion=2
		c_col=c_col*col_dens


	#same for H2
	if col_dens>1e23:
		h2_col=0.5
	else:
		h2_col=(10.0**random.uniform(-15.0,-0.3))

	param_dict["coldens"]=col_dens
	param_dict["ccol"]=c_col
	param_dict["h2col"]=h2_col*col_dens
	param_dict["rout"]=param_dict["coldens"]/param_dict["initialDens"]
	param_dict["ion"]=ion
	return param_dict

if __name__ == '__main__':
	models=[]

	av_factor=6.289E-22


	base_path=f"Dimensionality/raw/"
	df=pd.DataFrame()
	for i in range(1000):
		params=	{
			"columnFile":f"{base_path}{i}.dat",
			"outputFile":f"{base_path}{i}_full.dat",
			"finalTime":1.0e6,
			"outSpecies":12,
			"heatingFlag":True,
			"heatWriteFlag":False,
			"avFactor":av_factor}	

		params["initialTemp"]=random.uniform(10.0,1000.0)
		params["initialDens"]=10.0**random.uniform(2.0,6.0)
		params["radfield"]=10.0**random.uniform(0.0,5.0)
		params["zeta"]=10.0**random.uniform(0.0,5.0)
		params["metallicity"]=10.0**random.uniform(-3.0,0.0)

		params=generate_elements(params)
		params=generate_column_dens(params)

		models.append([i,params.copy()])
		params["model"]=i
		temp_df=pd.DataFrame(params,index=[1])
		df=df.append(temp_df,sort=False)
	df.to_csv(f"Dimensionality/models.csv",index=False)


	start=time.time()
	pool=Pool(24)
	print("mapping...")
	rel=pool.map(run_model,models)
	print(rel)
	pool.close()
	pool.join()
	end=time.time()
	print(start,end,end-start)
