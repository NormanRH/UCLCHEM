import uclchem
import numpy as np
import time
from multiprocessing import Pool
import pandas as pd

def run_model(a):
	i,cloud_size,h2_col,c_col,col_dens,model_type,rad_field=a
	params={
		"columnFile":f"Benchmarking/{model_type}/av_grid/{i:.0f}.dat",
		"initialTemp":100.0,
		"initialDens":1000.0,
		"finalTime":2.0e6,
		"radfield":rad_field,
		"outSpecies":12,
		"h2col":h2_col,
		"ccol":c_col,
		"coldens":col_dens,
		"rout":cloud_size,
		"fc":1.0e-4,
		"fhe":0.1,
		"fh":0.4,
		"fo":3.0e-4,
		"fr":0.0,
		"fmg":5.00e-06,
		"ion":2,
		"zeta":1.23}	

	outSpecies="H2O,H2,H,H+,HE+,C,C+,O,O+,CO,CO+,E-"
	uclchem.general(params,outSpecies,f"Benchmarking/{model_type}/av_grid/",f"{i:.0f}")
	return 1

if __name__ == '__main__':

	model_type="low_rad"
	if model_type=="low_rad":
		rad_field=2.89
	else:
		rad_field=850.0

	model_df=pd.read_csv(f"Benchmarking/{model_type}/model_tests.dat")
	print(model_df.head())
	models=[]
	with open(f"Benchmarking/{model_type}/av_grid.dat","w") as f:
		f.write("model,cloud_size,av,h2col,coldens\n")
		for i,row in model_df.iterrows():
			cloud_size=row["size"]/3.086e18
			h2_col=row["total_h_col"]
			c_col=row["total_c_col"]
			col_dens=row["total_col_dens"]

			# av conversion factor CHANGED TO MATCH UCL_PDR, should find most agreed value after benchmarking
			av=col_dens*1.6e-21

			print(cloud_size,av)
			models.append([i,cloud_size,h2_col,c_col,col_dens,model_type,rad_field])
			f.write(f"{i},{cloud_size},{av},{h2_col},{col_dens}\n")

	start=time.time()
	pool=Pool(24)
	print("mapping...")
	rel=pool.map(run_model,models)
	print(rel)
	pool.close()
	pool.join()
	end=time.time()
	print(start,end,end-start)
