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
speciesNameLists.append(("HS H2S #HS #H2S ").split())
speciesNameLists.append(("HCS OCS H2S2 #OCS #H2S2").split())
speciesNameLists.append(("NH NH2 #NH #NH2").split())
speciesNameLists.append(("NH3 HNC #NH3 #HNC").split())
speciesNameLists.append(("NO NO2 #NO #NO2").split())
speciesNameLists.append(("N HNCO #N #HNCO").split())
speciesNameLists.append(("HCN H2CN #HCN #H2CN").split())
speciesNameLists.append(("H2O HNO #H2O #HNO").split())
speciesNameLists.append(("H2O HNO #H2O #HNO").split())
speciesNameLists.append(("O O2 OH #O #O2 #OH").split())
speciesNameLists.append(("CH3 CH4 #CH3 #CH4").split())
speciesNameLists.append(("C3H2 CH3CCH #C3H2").split())
speciesNameLists.append(("CO CO2 #CO #CO2").split())
speciesNameLists.append(("SIC SIH2 #SIC #SIH4").split())

speciesiceNameLists = []
speciesiceNameLists.append(("SO SO2 S2").split())##SO #SO2
speciesiceNameLists.append(("HS H2S").split())# #HS #H2S
speciesiceNameLists.append(("HCS H2S2").split())# #HS #H2S
speciesiceNameLists.append(("NH NH2").split())##NH #NH2
speciesiceNameLists.append(("NH3 HNC").split())# #NH3 #HNC
speciesiceNameLists.append(("NO NO2").split())# #NO #NO2
speciesiceNameLists.append(("N HNCO").split())# #N #HNCO
speciesiceNameLists.append(("HCN H2CN").split())
speciesiceNameLists.append(("H2O HNO").split())
speciesiceNameLists.append(("O O2 OH").split())
speciesiceNameLists.append(("C CH CH2").split())# #C #CH #CH2
speciesiceNameLists.append(("CH3 CH4").split())# #CH3 #CH4
speciesiceNameLists.append(("C3H2 CH3CCH").split())
speciesiceNameLists.append(("CO CO2").split())# #CO #CO2
speciesiceNameLists.append(("SIC SIH2").split())# #SIC #HCN

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
        fig,axes=pf.plt.subplots(figsize=(16,9), sharey=True,num=p,clear=True)
        p=p+1
        #figcombo,axescombo=pf.plt.subplots(figsize=(16,9), sharey=True,num=p,clear=True)
        #p=p+1
#Separate out the phase 2 graph
        #figp2,axesp2=pf.plt.subplots(figsize=(16,9), sharey=True,num=p,clear=True)
        phase2startt=0
        phase2endt=0
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
            for l, s in enumerate(speciesNames):
                specfactor.append(s+"_" + str(factor))
      
      		#plot species and save to test.png, alternatively send dens instead of time.
            #axis,rtist0=pf.plot_species(speciesNames,time,abundances,axes,ls=Linestles[j],lab=lb)#ax=axes[i])
            #axis,rtist1=pf.plot_species(speciesNames,time,abundances,axes[1],ls=Linestles[j],lab=False)#ax=axes[i])
          
            #phase 2 data
            if factor == 1 : 
                collfile = columnpath + "phase2-fullCtrl.dat"
                #lb=True
                #title = "Ctrl"
            else: 
                collfile = columnpath + "phase2-full"+e[0]+ "-" + str(factor).replace('.','_')+".dat"
                #lb=False
                #title = "Varying " + e[0]#+str(factor)
      		#call read_uclchem. 
            time2,dens2,temp,abundances2=pf.read_uclchem(collfile,speciesNames)
          		#plot species returns the axis so we can further edit
            #stitch together the phase1 and phase2 series and plot
            timenp = np.asarray(time)
            timenp2 = np.asarray(time2)
            for l, s in enumerate(speciesNames):
                abundances[l].extend(abundances2[l])
            
            tnpp2 = timenp2 + timenp[-1]
            t = np.append(timenp,tnpp2)
            axis,rtist0=pf.plot_species(speciesNames,t,abundances,axes,ls=Linestles[j],lab=lb)#ax=axes[i])
            #axiscombo,rtist0=pf.plot_species(speciesNames,t,abundances,axescombo,ls=Linestles[j],lab=lb)#ax=axes[i])
            
            # plot phase2 on its own
            phase2startt=tnpp2[0]
            phase2endt=tnpp2[-1]
      		#plot species and save to test.png, alternatively send dens instead of time.
            #axisp2,rtist2=pf.plot_species(speciesNames,tnpp2,abundances2,axesp2,ls=Linestles[j],lab=lb)#ax=axes[i])
            
            
            
        axes.set(xscale='linear',yscale='log',ylim=(1e-18,1e-3),xlim=(1e0,timenp[-1]))
        axes.set_xlabel('Time / Myears')
        axes.ticklabel_format(axis='x',useMathText=True)
        axes.set_ylabel('X/H')
        axes.set_title(title+" phase1")
        axes.legend(loc='upper left',fontsize='small')
        fig.savefig("../VaryFromSolar/SepPlots"+str(switch)+"/phase1plot"+e[0]+str(m)+".png",dpi=300)
        
        axes.set(xscale='linear',yscale='log',ylim=(1e-18,1e-3),xlim=(1e0,6.5e6))#t[-1]))
        axes.set_title(title+" phase1&2")
        fig.savefig("../VaryFromSolar/SepPlots"+str(switch)+"/comboplot"+e[0]+str(m)+".png",dpi=300)

        axes.set(xscale='linear',yscale='log',ylim=(1e-18,1e-3),xlim=(phase2startt,phase2startt+0.3e6))
        axes.set_title(title+" phase2")
        fig.savefig("../VaryFromSolar/SepPlots"+str(switch)+"/phase2plot"+e[0]+str(m)+".png",dpi=300)

        
        #axescombo.set(xscale='linear',yscale='log',ylim=(1e-18,1e-3),xlim=(1e0,6.5e6))#t[-1]))
        #axescombo.set_xlabel('Time / Myears')
        #axescombo.ticklabel_format(axis='x',useMathText=True)
        #axescombo.set_ylabel('X/H')
        #axescombo.set_title(title+" phase1&2")
        #axescombo.legend(loc='upper left',fontsize='small')
        #figcombo.savefig("../VaryFromSolar/SepPlots"+str(switch)+"/comboplot"+e[0]+str(m)+".png",dpi=300)
        
        
        #axesp2.set(xscale='linear',yscale='log',ylim=(1e-18,1e-3),xlim=(phase2startt,phase2startt+0.3e6))
        #axisp2.set_ylabel('X/H')
        #axisp2.set_xlabel('Time / Myears')
        #axisp2.ticklabel_format(axis='x',useMathText=True)
        #axesp2.legend(loc='upper left',fontsize='small')
        #axesp2.set_title(title+" phase2")
        #figp2.savefig("../VaryFromSolar/SepPlots"+str(switch)+"/phase2plot"+e[0]+str(m)+".png",dpi=300)

        #axes[1].set(xscale='log',yscale='log',ylim=(1e-18,1e-3),xlim=(timenp[-1],2e7))
        #pf.plt.subplots_adjust(wspace=0.0)
        #the single plot per page version
        #axes.text(.02,0.98,e[0],horizontalalignment="left",verticalalignment="top",transform=axes.transAxes)
        #axes[3].text(.02,0.98,"Your Row",horizontalalignment="left",verticalalignment="top",transform=axes[3].transAxes)
        #pf.plt.clf()
        
#        p=p+1
#Separate out the phase 2 graph
#        fig,axes=pf.plt.subplots(figsize=(16,9), sharey=True,num=p,clear=True)
#        for j , factor in enumerate(varyfactor):
#            #pick species, any number is fine
#            #title = ""
#            lb=False
#            #phase 2 display
#            if factor == 1 : 
#                collfile = columnpath + "phase2-fullCtrl.dat"
#                lb=True
#                #title = "Ctrl"
#            else: 
#                collfile = columnpath + "phase2-full"+e[0]+ "-" + str(factor).replace('.','_')+".dat"
#                lb=False
#                #title = "Varying " + e[0]#+str(factor)
#      		#call read_uclchem. 
#            time2,dens2,temp,abundances2=pf.read_uclchem(collfile,speciesNames)
#            
#      		#plot species and save to test.png, alternatively send dens instead of time.
#            axis,rtist2=pf.plot_species(speciesNames,time2,abundances2,axes,ls=Linestles[j],lab=lb)#ax=axes[i])
#          
#          		#plot species returns the axis so we can further edit
#            #axis.set(xscale='log',ylim=(1e-18,1e-3),xlim=(1e0,6e6))
#            #axis.set_title(title)
#            #axis.legend(loc='best')
#            
#    
#        #axes[1].set_title(title+" phase1")
#        axes.set(xscale='log',yscale='log',ylim=(1e-18,1e-3),xlim=(1e0,2e6))
#        axis.set_ylabel('X/H')
#        axis.set_xlabel('Time / Myears')
#        axes.legend(loc='upper left',fontsize='small')
#        axes.set_title(title+" phase2")
#        #axes[1].legend(loc='best')
#        axes[0].text(.02,0.98,e[0],horizontalalignment="left",verticalalignment="top",transform=axes[0].transAxes)
#        pf.plt.subplots_adjust(wspace=0.0)
#        #the single plot per page version
#        #axes.text(.02,0.98,e[0],horizontalalignment="left",verticalalignment="top",transform=axes.transAxes)
#        #axes[3].text(.02,0.98,"Your Row",horizontalalignment="left",verticalalignment="top",transform=axes[3].transAxes)
#        fig.savefig("../VaryFromSolar/SepPlots"+str(switch)+"/phase2plot"+e[0]+str(m)+".png",dpi=300)
#        pf.plt.clf()

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
            axis.set(xscale='linear',yscale='log',ylim=(1e-18,1e-3),xlim=(1e0,3e6))
            axis.set_ylabel('X/H')
            axis.set_xlabel('Time / Myears')
            axis.ticklabel_format(axis='x',useMathText=True)
            axis.set_title(title+" static")
            axis.set_title(title)
            axis.legend(loc='upper left',fontsize='small')
            i=i+1
    
        #axes[0].text(.02,0.98,e[0],horizontalalignment="left",verticalalignment="top",transform=axes[0].transAxes)
        #the single plot per page version
#        axes.text(.02,0.98,e[0],horizontalalignment="left",verticalalignment="top",transform=axes.transAxes)
        
        #axes[3].text(.02,0.98,"Your Row",horizontalalignment="left",verticalalignment="top",transform=axes[3].transAxes)
        fig.savefig("../VaryFromSolar/SepPlots"+str(switch)+"/staticplot"+e[0]+str(m)+".png",dpi=300)
        pf.plt.clf()
