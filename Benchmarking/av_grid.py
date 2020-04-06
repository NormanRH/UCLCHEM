import uclchem
import numpy as np
import time
from multiprocessing import Pool


def run_model(a):
	i,cloud_size=a
	params={
		"columnFile":f"Benchmarking/high_rad/av_grid/{i:.0f}.dat",
		"initialTemp":100.0,
		"initialDens":1000.0,
		"radfield":850,
		"outSpecies":10,
		"rout":cloud_size
	}	

	outSpecies="H2O ,H2  ,H   ,CO  ,C+  ,O   ,CO+ ,CO2 ,#CO ,#H2O"
	uclchem.general(params,outSpecies,"Benchmarking/high_rad/av_grid/",f"{i:.0f}")


if __name__ == '__main__':
	models=[]
	with open("Benchmarking/high_rad/av_grid.dat","w") as f:
		f.write("model,cloud_size,av\n")
		for i,cloud_size in enumerate(np.logspace(-4.1,2.1,24)):
			av=3.086e18*cloud_size
			av=(av*1.0e3)/1.6e21
			print(cloud_size,av)
			models.append([i,cloud_size])
			f.write(f"{i},{cloud_size},{av}\n")

	start=time.time()
	pool=Pool(24)
	pool.map(run_model,models)
	pool.close()
	pool.join()
	end=time.time()
	print(start,end,end-start)
#uclchem.general(params,outSpecies)
