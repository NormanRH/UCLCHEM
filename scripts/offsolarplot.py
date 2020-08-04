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

speciesNames = ("SO SO2 S2 N NH NH2 NH3 HNC NO NO2 OCN HNCO HCS O2 H2O").split()

varyfactor = [0.25, 0.5, 1, 2, 4]
columnpath = "../VaryFromSolar/columnfiles/"

fig,axes=pf.plt.subplots(len(elements),len(varyfactor),figsize=(16,9))
axes=axes.flatten()
i=0

for k, e in enumerate(elements):
    iprev = i
    for j , factor in enumerate(varyfactor):
        #pick species, any number is fine
        title = ""
        if factor == 1 : 
            collfile = columnpath + "phase1-columnCtrl.dat"
            title = "Ctrl"
        else: 
            collfile = columnpath + "phase1-column"+e[0]+ "-" + str(factor).replace('.','_')+".dat"
            title = e[0]+str(factor)
  		#call read_uclchem. 
        time,dens,temp,abundances=pf.read_uclchem(collfile,speciesNames)
  
  
  		#plot species and save to test.png, alternatively send dens instead of time.
        axis=pf.plot_species(speciesNames,time,abundances,ax=axes[i])
      
      		#plot species returns the axis so we can further edit
        axis.set(xscale='log',ylim=(1e-15,1e-3),xlim=(1e0,6e6))
        axis.set_title(title)
        i=i+1

    axes[i-iprev].text(.02,0.98,e[0],horizontalalignment="left",verticalalignment="top",transform=axes[0].transAxes)
#axes[3].text(.02,0.98,"Your Row",horizontalalignment="left",verticalalignment="top",transform=axes[3].transAxes)
fig.savefig("../VaryFromSolar/Plots/speciesplot.png",dpi=300)
