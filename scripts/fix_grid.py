from __future__ import print_function
import uclchem
import numpy as np
import pandas as pd
import os
import time


models=pd.read_csv("output/grid.dat",delimiter=" ")

for i,row in models.iterrows():
    try:
        data=np.loadtxt(row["file"])
    except:
        print(row["file"]," failed")
        try:
            uclchem.temperaturetest(float(row["dens"]),float(row["temp"]),float(row["uv"]),1.0e6,row["file"])
        except:
            print("retry failed")