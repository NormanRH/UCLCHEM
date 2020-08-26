#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 20 21:10:09 2020

@author: norman
"""
#demonstration of plotfunctions. called from main UCLCHEM directory
#it reads full UCLCHEM output and saves a plot of the abudances of select species 

import numpy as np
import plotfunctions as pf
#Element ids and initial values
# This part can be substituted with any choice of grid
elements=[["C","fc",2.6e-04], ["N","fn",6.1e-05], ["O","fo",4.6e-04], ["S","fs",1.318e-05],
           ["Mg","fmg",3.981e-05], ["Si","fsi",1.0e-07],["Cl","fcl",3.162e-07],["F","ff",3.6e-08]]

#list the species we want graphing from the current data output of the model
#speciesNames = ("SO SO2 S2 N NH NH2 NH3 HNC NO NO2 OCN HNCO HCS O2 H2O").split()
speciesNameLists = []
speciesNameLists.append(("SO #SO SO2 #SO2 S2").split())
speciesNameLists.append(("HS H2S HCS #HS #H2S").split())
speciesNameLists.append(("NH NH2 #NH #NH2").split())
speciesNameLists.append(("NH3 HNC #NH3 #HNC").split())
speciesNameLists.append(("NO NO2 #NO #NO2").split())
speciesNameLists.append(("N HNCO #N #HNCO").split())
speciesNameLists.append(("C CH CH2 #C #CH #CH2").split())
speciesNameLists.append(("CH3 CH4 #CH3 #CH4").split())
speciesNameLists.append(("CO CO2 #CO #CO2").split())
speciesNameLists.append(("SIC HCN #SIC #HCN").split())

varyfactor = [0.25, 0.5, 1, 2, 4]
Linestles = [(0,(3,10,1,10)),(0,(3,5,1,5,1,5)),(0,()),(0,(3,5,1,5)),(0,(3,1,1,1))]

switch=1
columnpath = "../VaryFromSolar/outputfiles"+str(switch)+"/"

#fig,axes=pf.plt.subplots(len(elements),len(varyfactor),figsize=(16,9))
#fig,axes=pf.plt.subplots(figsize=(16,9))#len(elements),len(varyfactor),figsize=(16,9))
#axes=axes.flatten()
#i=0

for m , speciesNames in enumerate(speciesNameLists):
    p=0
    for k, e in enumerate(elements):
        i=0
        iprev = i
        specfactor = []
        title = "Varying " + e[0]
        #aiming to have 3 panes stacked vertically
        fig,axes=pf.plt.subplots(1,2,figsize=(16,9), sharey=True,num=p,clear=True)
        for j , factor in enumerate(varyfactor):
            #pick species, any number is fine
            #title = ""
            lb=False
            if factor == 1 : 
                collfile = columnpath + "phase1-fullCtrl.dat"
                #title = "Ctrl"
                lb=True
            else: 
                collfile = columnpath + "phase1-full"+e[0]+ "-" + str(factor).replace('.','_')+".dat"
                #title = "Varying " + e[0]#+str(factor)
                lb=False
      		#call read_uclchem. 
            time,dens,temp,abundances=pf.read_uclchem(collfile,speciesNames)
            for k, s in enumerate(speciesNames):
                specfactor.append(s+"_" + str(factor))
      
          
          		#plot species returns the axis so we can further edit
            #axis.set(xscale='log',ylim=(1e-18,1e-3),xlim=(1e0,6e6))
            #axis.set_title(title)
            #axis.legend(loc='best')
            
            #phase 2 display
            if factor == 1 : 
                collfile = columnpath + "phase2-fullCtrl.dat"
                #title = "Ctrl"
            else: 
                collfile = columnpath + "phase2-full"+e[0]+ "-" + str(factor).replace('.','_')+".dat"
                #title = "Varying " + e[0]#+str(factor)
      		#call read_uclchem. 
            time2,dens2,temp,abundances2=pf.read_uclchem(collfile,speciesNames)
            
            timenp = np.asarray(time)
            timenp2 = np.asarray(time2)
            for k, s in enumerate(speciesNames):
                abundances[k].extend(abundances2[k])
            
            t = np.append(timenp,(timenp2 + timenp[-1]))
            
            
      		#plot species and save to test.png, alternatively send dens instead of time.
            axis,rtist0=pf.plot_species(speciesNames,t,abundances,axes[0],ls=Linestles[j],lab=lb)#ax=axes[i])
            axis,rtist1=pf.plot_species(speciesNames,t,abundances,axes[1],ls=Linestles[j],lab=False)#ax=axes[i])
      		#plot species and save to test.png, alternatively send dens instead of time.
            #axis,rtist2=pf.plot_species(speciesNames,time2,abundances2,axes[2],ls=Linestles[j],lab=False)#ax=axes[i])
          
          		#plot species returns the axis so we can further edit
            #axis.set(xscale='log',ylim=(1e-18,1e-3),xlim=(1e0,6e6))
            #axis.set_title(title)
            #axis.legend(loc='best')
            
    
        axes[0].set(xscale='log',yscale='log',ylim=(1e-18,1e-3),xlim=(1e0,timenp[-1]))
        axes[0].set_xlabel('Time / years')
        axes[0].set_ylabel('X/H')
        axes[0].set_title(title+" phase1")
        axes[0].legend(loc='upper left',fontsize='small')
        axes[1].set(xscale='log',yscale='log',ylim=(1e-18,1e-3),xlim=(timenp[-1],t[-1]))
        axes[1].set_title("Zoom in on phase 2")
        #axes[1].set_title(title+" phase1")
        #axes[2].set(xscale='log',yscale='log',ylim=(1e-18,1e-3),xlim=(1e2,6e6))
        #axes[2].set_title(title+" phase2")
        #axes[1].legend(loc='best')
#        axes[0].text(.02,0.98,e[0],horizontalalignment="left",verticalalignment="top",transform=axes[0].transAxes)
        pf.plt.subplots_adjust(wspace=0.0)
        #the single plot per page version
        #axes.text(.02,0.98,e[0],horizontalalignment="left",verticalalignment="top",transform=axes.transAxes)
        #axes[3].text(.02,0.98,"Your Row",horizontalalignment="left",verticalalignment="top",transform=axes[3].transAxes)
        fig.savefig("../VaryFromSolar/Plots"+str(switch)+"/speciesplot"+e[0]+str(m)+".png",dpi=300)
        #pf.plt.clf()
    p=p+1
    for k, e in enumerate(elements):
        i=0
        iprev = i
        specfactor = []
        #aiming to have 3 panes stacked vertically
        fig,axes=pf.plt.subplots(figsize=(16,9),num=p,clear=True)
        for j , factor in enumerate(varyfactor):
            #pick species, any number is fine
            title = ""
            lb=False
            
            #static model plot
            if factor == 1 : 
                collfile = columnpath + "static-fullCtrl.dat"
                title = "Ctrl"
                lb = True
            else: 
                collfile = columnpath + "static-full"+e[0]+ "-" + str(factor).replace('.','_')+".dat"
                title = "Varying " + e[0]#+str(factor)
                lb = False
      		#call read_uclchem. 
            time3,dens3,temp3,abundances3=pf.read_uclchem(collfile,speciesNames)
      
      		#plot species and save to test.png, alternatively send dens instead of time.
            axis,rtist=pf.plot_species(speciesNames,time3,abundances3,axes,ls=Linestles[j],lab=lb)#ax=axes[i])
          
          		#plot species returns the axis so we can further edit
            axis.set(xscale='log',yscale='log',ylim=(1e-18,1e-3),xlim=(1e0,time3[-1]))
            axis.set_ylabel('X/H')
            axis.set_xlabel('Time / years')
            axis.set_title(title+" static")
            axis.set_title(title)
            axis.legend(loc='upper left',fontsize='small')
            i=i+1
    
        #axes[0].text(.02,0.98,e[0],horizontalalignment="left",verticalalignment="top",transform=axes[0].transAxes)
        #the single plot per page version
#        axes.text(.02,0.98,e[0],horizontalalignment="left",verticalalignment="top",transform=axes.transAxes)
        
        #axes[3].text(.02,0.98,"Your Row",horizontalalignment="left",verticalalignment="top",transform=axes[3].transAxes)
        fig.savefig("../VaryFromSolar/Plots"+str(switch)+"/staticplot"+e[0]+str(m)+".png",dpi=300)
        #pf.plt.clf()
