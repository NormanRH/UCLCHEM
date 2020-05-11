from plotting_functions import *

base_path = 'Benchmarking/'
for model_type in ["10_1e5.5","1e5_1e5.5","10_1e3","1e5_1e3","low_rad_fixed","high_cr"]:#,"fixed_cooling"]:
	path=f"{base_path}raw/{model_type}"
	uclchem_df=create_uclchem_df(path)
	uclchem_df.to_csv(f"{base_path}/uclchem/{model_type}.csv")
