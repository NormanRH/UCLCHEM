import uclchem
import numpy as np
import time
from multiprocessing import Pool
import pandas as pd
from os import path,mkdir

def run_model(a):
	i,cloud_size,h2_col,c_col,col_dens,model_type,rad_field,init_temp,dens,heating_flag,av_factor,zeta,cooling,heating=a
	params={
		"columnFile":f"Benchmarking/{model_type}/Equilibrium/av_grid/{i:.0f}.dat",
		"initialTemp":init_temp,
		"initialDens":dens,
		"finalTime":2.0e6,
		"radfield":rad_field,
		"outSpecies":12,
		"h2col":h2_col,
		"ccol":c_col,
		"coldens":col_dens,
		"rout":cloud_size,
		"fc":1.0000000e-4,
		"fhe":0.1,
		"fh":0.4,
		"fo":3.0e-4,
		"fr":0.0,
		"fmg":5.00e-06,
		"ion":2,
		"zeta":zeta,
		"heatingFlag":heating_flag,
		"avFactor":av_factor,
		"cooling":cooling,
		"heating":heating}	
	print(i)
	outSpecies="H2O,H2,H,H+,HE+,C,C+,O,O+,CO,CO+,E-"
	uclchem.equilibrium(params,outSpecies,f"Benchmarking/{model_type}/Equilibrium/av_grid/",f"{i:.0f}")
	return 1

if __name__ == '__main__':
	models=[]

	av_converts={"10_1e3":6.289E-22,
				"10_1e5.5":6.289E-22,
				"1e5_1e3":6.289E-22,
				"1e5_1e5.5":6.289E-22,
				"fixed_cooling":6.289E-22,
				"fixed_heating":6.289E-22,
				"low_rad":1.6e-21,
				"low_rad_fixed":1.6e-21,
				"high_cr":1.6e-21}

	zetas={"10_1e3":3.84,
		"10_1e5.5":3.84,
		"1e5_1e3":3.84,
		"1e5_1e5.5":3.84,
		"fixed_cooling":3.84,
		"fixed_heating":3.84,
		"low_rad":1.23,
		"low_rad_fixed":1.23,
		"high_cr":123.0
	}


	for model_type in ["fixed_heating"]:#"fixed_cooling","10_1e5.5","1e5_1e5.5","10_1e3","1e5_1e3","low_rad_fixed","high_cr"
		fixed_temp=("low_rad_fixed" in model_type)
		fixed_cooling=("fixed_cooling" in model_type)
		fixed_heating=("fixed_heating" in model_type)

		base_path=f"Benchmarking/{model_type}/Equilibrium/"
		if not path.exists(base_path):
			mkdir(base_path)

		if not path.exists(base_path+"av_grid/"):
			mkdir(base_path+"av_grid/")
		if not path.exists(base_path+"av_grid/cooling/"):
			mkdir(base_path+"av_grid/cooling/")
		if not path.exists(base_path+"av_grid/heating/"):
			mkdir(base_path+"av_grid/heating/")


		#prop.out has radiation field for each particle. first one is the surface field we need
		model_df=pd.read_csv(f"{base_path}{model_type}.prop.out",nrows=5)
		#but we need it in Habing for UCLCHEM
		rad_field=model_df.loc[1,"chi"]*1.71

		model_df=pd.read_csv(f"{base_path}model_tests.dat")
		model_df=model_df.merge(pd.read_csv(f"{base_path}{model_type}.cool.out",usecols=["Particle","Total"],delimiter=" "),
								on="Particle",how="left")
		model_df=model_df.rename({"Total":"totalCooling"},axis=1)
		model_df=model_df.merge(pd.read_csv(f"{base_path}{model_type}.heat.out",usecols=["Particle","Total"],delimiter=" "),
								on="Particle",how="left")
		with open(f"{base_path}av_grid.dat","w") as f:
			f.write("model,cloud_size,av,h2col,coldens\n")
			for i,row in model_df.iterrows():
				cloud_size=row["size"]/3.086e18
				h2_col=row["total_h_col"]
				c_col=row["total_c_col"]
				col_dens=row["total_col_dens"]
				dens=row["n_H"]

				if fixed_temp:
					initialTemp=row["T_gas"]
					heatingFlag=False
				
				else:
					initialTemp=row["T_gas"]
					heatingFlag=True


				if fixed_cooling:
					cooling=row["totalCooling"]
				else:
					cooling=0.0

				if fixed_heating:
					heating=row["Total"]
				else:
					heating=0.0

				# av conversion factor CHANGED TO MATCH UCL_PDR, should find most agreed value after benchmarking
				av=col_dens*av_converts[model_type]
				zeta=zetas[model_type]
				models.append([i,cloud_size,h2_col,c_col,col_dens,model_type,rad_field,initialTemp,dens,heatingFlag,av_converts[model_type],zeta,cooling,heating])
				f.write(f"{i},{cloud_size},{av},{h2_col},{col_dens}\n")

	start=time.time()
	pool=Pool(4)
	print("mapping...")
	rel=pool.map(run_model,models)
	print(rel)
	pool.close()
	pool.join()
	end=time.time()
	print(start,end,end-start)
