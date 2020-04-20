import uclchem
import numpy as np
import time
from multiprocessing import Pool
import pandas as pd

def run_model(a):
	i,cloud_size,h2_col,c_col,col_dens,model_type,rad_field,init_temp,heating_flag,av_factor,zeta=a	
	params={
		#"columnFile":f"Benchmarking/{model_type}/av_grid/{i:.0f}.dat",
		"columnFile":"test/test.dat",
		"initialTemp":init_temp,
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
		"zeta":zeta,
		"heatingFlag":heating_flag,
		"avFactor":av_factor}	

	outSpecies="H2O,H2,H,H+,HE+,C,C+,O,O+,CO,CO+,E-"
	uclchem.general(params,outSpecies,"test/","1")
	return 1


av_converts={"10_1e3":6.289E-22,
			"10_1e5.5":6.289E-22,
			"1e5_1e3":6.289E-22,
			"1e5_1e5.5":6.289E-22,
			"low_rad":1.6e-21,
			"low_rad_fixed":1.6e-21,
			"high_cr":1.6e-21}

zetas={"10_1e3":3.84,
	"10_1e5.5":3.84,
	"1e5_1e3":3.84,
	"1e5_1e5.5":3.84,
	"low_rad":1.23,
	"low_rad_fixed":1.23,
	"high_cr":123.0
}


model_type="1e5_1e5.5"
fixed="fixed" in model_type
base_path=f"Benchmarking/{model_type}/"

#prop.out has radiation field for each particle. first one is the surface field we need
model_df=pd.read_csv(f"{base_path}{model_type}.prop.out",nrows=5)
#but we need it in Habing for UCLCHEM
rad_field=model_df.loc[1,"chi"]*1.7

i=3
model_df=pd.read_csv(f"Benchmarking/{model_type}/model_tests.dat")
row=model_df.iloc[i]
cloud_size=row["size"]/3.086e18
h2_col=row["total_h_col"]
c_col=row["total_c_col"]
col_dens=row["total_col_dens"]

if fixed:
	initialTemp=row["T_gas"]
	heatingFlag=False
else:
	initialTemp=50.0
	heatingFlag=True

# av conversion factor CHANGED TO MATCH UCL_PDR, should find most agreed value after benchmarking
av=col_dens*av_converts[model_type]
zeta=zetas[model_type]

success=run_model([i,cloud_size,h2_col,c_col,col_dens,model_type,rad_field,initialTemp,heatingFlag,av_converts[model_type],zeta])
print(success)