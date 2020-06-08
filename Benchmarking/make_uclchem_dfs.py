from plotting_functions import *

base_path = 'Benchmarking/'
for model_type in ["low_rad"]:#,"fixed_cooling"]:"low_metallicity","10_linear_increase","10_linear_decrease","10_1e5.5","1e5_1e5.5","10_1e3","1e5_1e3","high_cr",
	path=f"Benchmarking/reduced-raw/{model_type}"
	uclchem_df=create_uclchem_df(path)
	uclchem_df.to_csv(f"{base_path}/uclchem-reduced/{model_type}.csv")
