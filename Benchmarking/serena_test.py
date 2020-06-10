import uclchem
import numpy as np
import time
from multiprocessing import Pool
import pandas as pd
from os import path,mkdir

params={
	"columnFile":f"test/sv_col.dat",
	"outputFile":f"test/sv_full.dat",
	"initialTemp":100.0,
	"initialDens":1.0e5,
	"finalTime":2.0e6,
	"radfield":100.0,
	"outSpecies":12,
	"h2col":3.1800652239385913e+22,
	"ccol":4.0071298998096024e+16,
	"coldens":6.3603163110564155e+22,
	"rout":2.01130713e+17,
	 "fc":1.0e-4,
	 "fhe":0.1,
	 "fh":0.4,
	 "fo":3.0e-4,
	 "fn":0.0,
	 "fs":0.0,#1.318e-8,
	 "fcl":0.0,#3.162e-10,
	 "fsi":0.0,#1.0e-10,
	 "fmg":5.00e-06,
	"ion":2,
	"zeta":10000.0,
	"fr":0.0,
	"metallicity":1.0,
	"heatingFlag":True,
	"heatWriteFlag":False,
	"avFactor":6.289E-22}

start=time.time()
outSpecies="H2O,H2,H,H+,HE+,C,C+,O,O+,CO,CO+,E-"
uclchem.general(params,outSpecies,f"test/","1")
end=time.time()

print(f"Completed in {end-start} seconds")