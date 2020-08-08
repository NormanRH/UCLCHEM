#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 20 21:10:09 2020

@author: norman
"""
#demonstration of plotfunctions. called from main UCLCHEM directory
#it reads full UCLCHEM output and saves a plot of the abudances of select species 

import plotfunctions as pf
#Element ids and initial values
# This part can be substituted with any choice of grid
elements=[["C","fc",2.6e-04], ["N","fn",6.1e-05], ["O","fo",4.6e-04], ["S","fs",1.318e-05],
           ["Mg","fmg",3.981e-05], ["Si","fsi",1.0e-07],["Cl","fcl",3.162e-07],["F","ff",3.6e-08]]

#list the species we want graphing from the current data output of the model
#speciesNames = ("SO SO2 S2 N NH NH2 NH3 HNC NO NO2 OCN HNCO HCS O2 H2O").split()
speciesNames = ("SO SO2 S2").split()

varyfactor = [0.25, 0.5, 1, 2, 4]
Linestles = [(0,(3,10,1,10)),(0,(3,5,1,5,1,5)),(0,()),(0,(3,5,1,5)),(0,(3,1,1,1))]
columnpath = "../VaryFromSolar/outputfiles/"

#fig,axes=pf.plt.subplots(len(elements),len(varyfactor),figsize=(16,9))
#fig,axes=pf.plt.subplots(figsize=(16,9))#len(elements),len(varyfactor),figsize=(16,9))
#axes=axes.flatten()
#i=0

for k, e in enumerate(elements):
    i=0
    iprev = i
    specfactor = []
    fig,axes=pf.plt.subplots(figsize=(16,9),num=k,clear=True)
    for j , factor in enumerate(varyfactor):
        #pick species, any number is fine
        title = ""
        if factor == 1 : 
            collfile = columnpath + "phase1-fullCtrl.dat"
            title = "Ctrl"
        else: 
            collfile = columnpath + "phase1-full"+e[0]+ "-" + str(factor).replace('.','_')+".dat"
            title = "Varying " + e[0]#+str(factor)
  		#call read_uclchem. 
        time,dens,temp,abundances=pf.read_uclchem(collfile,speciesNames)
        for k, s in enumerate(speciesNames):
            specfactor.append(s+"_strfactor")
  
  		#plot species and save to test.png, alternatively send dens instead of time.
        axis=pf.plot_species(speciesNames,time,abundances,axes,ls=Linestles[j])#ax=axes[i])
      
      		#plot species returns the axis so we can further edit
        axis.set(xscale='log',ylim=(1e-15,1e-3),xlim=(1e0,6e6))
        axis.set_title(title)
        i=i+1

    #axes[i-iprev].text(.02,0.98,e[0],horizontalalignment="left",verticalalignment="top",transform=axes[0].transAxes)
    #the single plot per page version
    axes.text(.02,0.98,e[0],horizontalalignment="left",verticalalignment="top",transform=axes.transAxes)
    #axes[3].text(.02,0.98,"Your Row",horizontalalignment="left",verticalalignment="top",transform=axes[3].transAxes)
    fig.savefig("../VaryFromSolar/Plots/speciesplot"+e[0]+".png",dpi=300)
    #fig.clear()
