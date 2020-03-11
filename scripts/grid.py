from __future__ import print_function
import uclchem
import numpy as np
from multiprocessing import Pool
import os
import time

#uclchem(final_dens,max_temp,shock_vel,phase_flag,outFile,startFile)
def runuclchem(a):
    modelNo,initial_dens,uv,fin_time,temp,outFile=a
    uclchem.temperaturetest(initial_dens,temp,uv,fin_time,outFile)

print(__name__)

#Set up grid by writing all parameter combinations to a file
#and creating a list of lists of parameters
#pass all lists to pool.map() to be done in parallel
#remember to do phase 1s first if doing multiple phases
if __name__ == '__main__':
    outputFolder="output/"
    f=open(outputFolder+"grid.dat","w")
    f.write("model dens temp uv file")
    densities = [100,1.0e4,1.0e6]
    temperatures=[10.0,30.0,50.0,100.0,300.0,500.0,1000.0]#np.logspace(0.5,10)
    uvs=[1.0,10.0,100.0]
    modelNo=1

    models=[]
    for init_dens in densities:
        print(init_dens)
        for temp in temperatures:
            for uv in uvs:
                outFile=outputFolder+"data/{0:.0f}".format(modelNo)
                if not os.path.isfile(outFile):
                    models.append([modelNo,init_dens,uv,1e6,temp,outFile])
                f.write("{0} {1} {2} {3} {4}\n".format(modelNo,init_dens,temp,uv,outFile))
                modelNo+=1

    f.close()
    start=time.time()
    pool=Pool(8)
    pool.map(runuclchem,models)
    pool.close()
    pool.join()
    end=time.time()
    with open("grid_time.dat","w") as f:
        f.write("grid in {0:.3f} seconds".format(end-start))