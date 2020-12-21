#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 18 11:30:03 2020

@author: norman hansen
"""

#Acknowledgement to Marcus Keil and Jon Holdship 13/03/2020
#Examples of a simple grid of models run in parallel

from __future__ import print_function
import os
import uclchem
import numpy as np
import pandas as pd
from multiprocessing import Pool, Process
import time

#uclchem general takes a dictionary of parameters where outSpecies=number of outspecies
#and a string outspeciesin which is a delimited list of species
#this wrapper helps Pool to work and also converts a delimited list in the dictionary into a number
#and send list separately as required.
def run_uclchem(param_dict):
    outSpecies = (param_dict['outSpecies'])
    param_dict['outSpecies'] = len(outSpecies.split())
    
    #this will run uclchem to final time step and return an array of abundances
    #arrays cannot be variably sized with f2py so you need to set length of output array in src/wrap.f90
    #abunds=uclchem.wrap.run_model_for_abundances(dictionary=param_dict, outspeciesin=outSpecies)
    abunds=uclchem.wrap.run_model_to_file(dictionary=param_dict, outspeciesin=outSpecies)
    
    #altenatively, you can run uclchem to file as you normally would
    #just put columnfile and/or output file in the dictionary
    #uclchem.wrap.run_model_to_file(dictionary, outSpeciesIn)
    return abunds

#path for output direcoties
model = "UCLCHEM"
switch=1
basepath="../VaryFromSolar"+model+"/"
if os.path.isdir(basepath) is not True:
    os.mkdir(basepath)

intermediatepath="../VaryFromSolar"+model+"/intermediatefiles"+str(switch)+"/"
if os.path.isdir(intermediatepath) is not True:
    os.mkdir(intermediatepath)
outputpath = "../VaryFromSolar"+model+"/outputfiles"+str(switch)+"/"
if os.path.isdir(outputpath) is not True:
    os.mkdir(outputpath)
columnpath = "../VaryFromSolar"+model+"/columnfiles"+str(switch)+"/"
if os.path.isdir(columnpath) is not True:
    os.mkdir(columnpath)
#basic set of parameters we'll use for this grid. 
#problem of choice of species and keeping track, need representative set of species but can we get 
#UCLChem to write files for susequent analysis we really want to output the 
#Start with this set as a static model and ignore a phase 2
ParameterDictP1 = {    "phase": 1, 
                       "switch": switch, 
                       "collapse": 1, 
                       "readAbunds": 0, 
                       "writeStep": 1,#"fn": 6.1e-9,
                       "outSpecies": 'SO SO2 S2 N NH NH2 NH3 HNC NO NO2 OCN HNCO HCS O2 H2O',
                       "fc":2.6e-04,
                       "fn":6.1e-05,
                       "fo":4.6e-04,
                       "fs":1.318e-05,
                       "fmg":3.981e-05,
                       "fsi":1.0e-07,
                       "fcl":3.162e-07,
                       "ff":3.6e-08,
                       "initialTemp": 10.0,
                       "initialDens": 1e2,
                       "finalDens":1.0e5,
                       "finalTime":1.0e7,#default final time
                       "baseAv":10}
#Phase2 after the above we continue the chemistry of collapse
ParameterDictP2 = {    "phase": 2, 
                       "switch": 0,#switch, 
                       "collapse": 0, 
                       "readAbunds": 1, #"writeStep": 1,
                       "tempindx": 3,
                       "fr": 0.0,
                       "outSpecies": 'SO SO2 S2 N NH NH2 NH3 HNC NO NO2 OCN HNCO HCS O2 H2O',
                       "fc":2.6e-04,
                       "fn":6.1e-05,
                       "fo":4.6e-04,
                       "fs":1.318e-05,
                       "fmg":3.981e-05,
                       "fsi":1.0e-07,
                       "fcl":3.162e-07,
                       "ff":3.6e-08,
                       "initialTemp": 10.0,
                       "initialDens": 1.0e5,
                       "finalTime":2.0e6,
                       "baseAv":10}

#Static model
ParameterDictP3 = {    "phase": 1, 
                       "switch": 0,#switch,#final Dens here is the same as start so need to run for time 
                       "collapse": 0, 
                       "readAbunds": 0, 
                       "writeStep": 1,
                       "tempindx": 1,
                       "fr": 0.0,
                       "outSpecies": 'SO SO2 S2 N NH NH2 NH3 HNC NO NO2 OCN HNCO HCS O2 H2O',
                       "fc":2.6e-04,
                       "fn":6.1e-05,
                       "fo":4.6e-04,
                       "fs":1.318e-05,
                       "fmg":3.981e-05,
                       "fsi":1.0e-07,
                       "fcl":3.162e-07,
                       "ff":3.6e-08,
                       "initialTemp": 10.0,
                       "initialDens": 1.0e4,
                       "finalDens": 1.0e4,
                       "finalTime":1.0e8,
                       "baseAv":10}

#Element ids and initial values
# This part can be substituted with any choice of grid

ModelelementsDict = { "UCLCHEM":[["C","fc",2.6e-04],
                          ["N","fn",6.1e-05],
                          ["O","fo",4.6e-04],
                          ["S","fs",1.318e-05],
                          ["Mg","fmg",3.981e-05],
                          ["Si","fsi",1.0e-07],
                          ["Cl","fcl",3.162e-07],
                          ["F","ff",3.6e-08]],
            
            "Wakelam":[["C","fc",1.2e-04],
                      ["N","fn",7.6e-05],
                      ["O","fo",2.56e-04],
                      ["S","fs",8.0e-08],
                      ["Mg","fmg",0.0],
                      ["Si","fsi",0.0],
                      ["Cl","fcl",1.8e-08],
                      ["F","ff",6.68e-09]],
            
            
            "NavAlm":[["C","fc",1.7e-04],
                      ["N","fn",6.2e-05],
                      ["O","fo",2.4e-04],
                      ["S","fs",1.5e-05],
                      ["Mg","fmg",7.0e-09],
                      ["Si","fsi",8.0e-09],
                      ["Cl","fcl",1.0e-09],
                      ["F","ff",6.68e-09]]
            }
elements = ModelelementsDict[model]
varyfactor = [0.25, 0.5, 1, 2, 4]#1 is the control
# here we just combine various initial and final densities into an easily iterable array
#initialTempGrid = np.linspace(50, 300, 6)
#crgrid = np.logspace(1, 5, 5)
#parameterSpace = np.asarray(np.meshgrid(initialTempGrid, crgrid)).reshape(2, -1)
v=open("varyfactor.txt","w")
print(varyfactor,v)
v.close()

e=pd.DataFrame(elements)
e.to_csv("elements.csv",index=False)
#then we loop through parameters, update parameter dictionary and run
#however, to use Pool we just store dictionaries for now
models=[]
models2=[]
models3=[]
for k, e in enumerate(elements):
    ParameterDictP1[e[1]]=e[2]
    ParameterDictP2[e[1]]=e[2]
    ParameterDictP3[e[1]]=e[2]
#we'll number our models so store the parameters we used in a list
for k, e in enumerate(elements):
    for j , factor in enumerate(varyfactor):
        if factor == 1:
            continue #1 is a common nothing changed so do once at the end
        paramDict1=ParameterDictP1.copy()
        paramDict2=ParameterDictP2.copy()
        paramDict3=ParameterDictP3.copy()
        paramDict1[e[1]] = e[2] * factor #only change the element we are running for now
        abundfile= intermediatepath + "startcollapse"+e[0]+ "-" + str(factor).replace('.','_') +".dat"
        paramDict1["abundFile"] = abundfile
        paramDict2["abundFile"] = abundfile #feed into phase 2 from phase1
        paramDict3[e[1]] = e[2] * factor
        paramDict3["abundFile"] = intermediatepath + "staticabund"+e[0]+ "-" + str(factor).replace('.','_') +".dat"
        paramDict1["outputFile"]= outputpath + "phase1-full"+e[0]+ "-" + str(factor).replace('.','_') +".dat"
        paramDict2["outputFile"]= outputpath + "phase2-full"+e[0] + "-" +str(factor).replace('.','_') +".dat"
        paramDict3["outputFile"]= outputpath + "static-full"+e[0] + "-" +str(factor).replace('.','_') +".dat"
        paramDict1["columnFile"]= columnpath + "phase1-column"+e[0]+ "-" + str(factor).replace('.','_')+".dat"
        paramDict2["columnFile"]= columnpath + "phase2-column"+e[0] + "-" +str(factor).replace('.','_')+".dat"
        paramDict3["columnFile"]= columnpath + "static-column"+e[0] + "-" +str(factor).replace('.','_')+".dat"
        #paramDict["zeta"] = parameterSpace[1,k]
        models.append(paramDict1)
        models2.append(paramDict2)
        models3.append(paramDict3)

#Control run with no deviation
paramDict1=ParameterDictP1.copy()
paramDict2=ParameterDictP2.copy()
paramDict3=ParameterDictP3.copy()
abundfile= intermediatepath + "startcollapseCtrl.dat"
paramDict1["abundFile"] = abundfile
paramDict2["abundFile"] = abundfile #feed into phase 2 from phase1
paramDict3["abundFile"] = intermediatepath + "staticCtrl.dat"
paramDict1["outputFile"]= outputpath + "phase1-fullCtrl.dat"
paramDict2["outputFile"]= outputpath + "phase2-fullCtrl.dat"
paramDict3["outputFile"]= outputpath + "static-fullCtrl.dat"
paramDict1["columnFile"]= columnpath + "phase1-columnCtrl.dat"
paramDict2["columnFile"]= columnpath + "phase2-columnCtrl.dat"
paramDict3["columnFile"]= columnpath + "static-columnCtrl.dat"
#paramDict["zeta"] = parameterSpace[1,k]
models.append(paramDict1)
models2.append(paramDict2)
models3.append(paramDict3)

m=pd.DataFrame(models)
m.to_csv("models.csv",index=False)
#use pool.map to run each dictionary throuh our helper function
start=time.time()
pool=Pool(16)
pool.map(run_uclchem,models)
#result=pool.map(run_uclchem,models)
#result=np.asarray(result)
pool.close()
pool.join()

pool = Pool(16)
pool.map(run_uclchem,models2)

#result2=pool.map(run_uclchem,models2)
#result2 = np.asarray(result2)
pool.close()
pool.join()

#static model
pool = Pool(16)
pool.map(run_uclchem,models3)

#result2=pool.map(run_uclchem,models2)
#result2 = np.asarray(result2)
pool.close()
pool.join()

#df=pd.DataFrame({#"Element":elements[0,:],"Mult":varyfactor[:],
#                "SO":result[:,0],"H3O+":result[:,1],"NH3":result[:,2],"HNC":result[:,3]})

#df.to_csv("grid.csv",index=False)
end=time.time()
end=(end-start)/60.0
print(f"grid in {end} minutes")

