from plotting_functions import *

base_path = 'Benchmarking/'
for model_type in ["low_metallicity"]:#,"fixed_cooling"]:
	path=f"Benchmarking/raw/{model_type}"
	uclchem_df=create_uclchem_df(path)
	uclchem_df.to_csv(f"{base_path}/uclchem/{model_type}.csv")
