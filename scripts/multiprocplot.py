#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 20 21:10:09 2020

@author: norman
"""
#demonstration of plotfunctions. called from main UCLCHEM directory
#it reads full UCLCHEM output and saves a plot of the abudances of select species 
import os
import multiprocessing as mp
import numpy as np
import plotfunctions as pf
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.lines as mlines
#Element ids and initial values
# This part can be substituted with any choice of grid
elements=[["C","fc",2.6e-04], ["N","fn",6.1e-05], ["O","fo",4.6e-04], ["S","fs",1.318e-05],
           ["Mg","fmg",3.981e-05], ["Si","fsi",1.0e-07],["Cl","fcl",3.162e-07],["F","ff",3.6e-08]]

imgsize = {"A7":[(3.5,2),3.0,4.0,4.0,0.5,4],
           "A6":[(5.2,3.5),'xx-small','x-small','small',1.0,7],
           "A5":[(6.8,4.5),'x-small','small','medium',1.5,8],
           "A4":[(10,6.8),'small','medium','large',2.0,9]}

#list the species we want graphing from the current data output of the model
#speciesNames = ("SO SO2 S2 N NH NH2 NH3 HNC NO NO2 OCN HNCO HCS O2 H2O").split()
#
speciesNameLists = []
speciesNameLists.append(("SO SO2 S2 #SO2 #SO").split())
speciesNameLists.append(("HS H2S #HS #H2S").split())
speciesNameLists.append(("HCS OCS H2S2 #OCS #H2S2").split())
speciesNameLists.append(("NH NH2 #NH #NH2").split())
speciesNameLists.append(("NH3 HNC #NH3 #HNC").split())
speciesNameLists.append(("NO NO2 #NO #NO2").split())
speciesNameLists.append(("N HNCO #N #HNCO").split())
speciesNameLists.append(("HCN H2CN #HCN #H2CN").split())
speciesNameLists.append(("H2O HNO #H2O #HNO").split())
speciesNameLists.append(("O O2 OH #O #O2 #OH").split())
speciesNameLists.append(("CH C CH2 #C #CH #CH2").split())# #C #CH #CH2
speciesNameLists.append(("CH3 CH4 #CH3 #CH4").split())
speciesNameLists.append(("C3H2 CH3CCH #C3H2").split())
speciesNameLists.append(("CO CO2 #CO #CO2").split())
speciesNameLists.append(("SIC SIH2 #SIC #SIH4").split())
speciesNameLists.append(("CL HCL MG #HCL").split())
speciesNameLists.append(("SIS SIC4 H2SIO SIO").split())
speciesNameLists.append(("MG MG+ E-").split())
speciesNameLists.append(("S S+ E-").split())
speciesNameLists.append(("C C+ E-").split())

speciesNoiceNameLists = []
speciesNoiceNameLists.append(("SO SO2 S2").split())##SO #SO2
speciesNoiceNameLists.append(("HS H2S").split())# #HS #H2S
speciesNoiceNameLists.append(("HCS OCS").split())# #HS #H2S"HCS OCS H2S2"
speciesNoiceNameLists.append(("NH NH2").split())##NH #NH2
speciesNoiceNameLists.append(("NH3 HNC").split())# #NH3 #HNC
speciesNoiceNameLists.append(("NO NO2").split())# #NO #NO2
speciesNoiceNameLists.append(("N HNCO").split())# #N #HNCO
speciesNoiceNameLists.append(("HCN H2CN").split())
speciesNoiceNameLists.append(("H2O HNO").split())
speciesNoiceNameLists.append(("O O2 OH").split())
speciesNoiceNameLists.append(("CH C CH2").split())# #C #CH #CH2
speciesNoiceNameLists.append(("CH3 CH4").split())# #CH3 #CH4
speciesNoiceNameLists.append(("C3H2 CH3CCH").split())
speciesNoiceNameLists.append(("CO CO2").split())# #CO #CO2
speciesNoiceNameLists.append(("SIC SIH2").split())# #SIC #HCN
speciesNoiceNameLists.append(("CL HCL MG").split())# #SIC #HCN
speciesNoiceNameLists.append(("SIS SIC4 H2SIO SIO").split())
speciesNoiceNameLists.append(("MG MG+ E-").split())
speciesNoiceNameLists.append(("S S+ E-").split())
speciesNoiceNameLists.append(("C C+ E-").split())

varyfactor = [0.25, 0.5, 1, 2, 4]
varyfactorstr = ["0.25", "0.5", "1", "2", "4"]
#Linestles = [(0,(4,2,1,2,1,2,1,2,1,2)),(0,(4,2,1,2,1,2)),(0,()),(0,(1,2,4,2,4,2)),(0,(1,2,4,2,4,2,4,2,4,2))]
#Linestles = [(0,(3,10,1,10)),(0,(3,5,1,5,1,5)),(0,()),(0,(3,5,1,5)),(0,(3,1,1,1))]
Linestles = [(0,(1,2)),(0,(1,1)),(0,()),(0,(2,1,1,1)),(0,(2,1))]

bulk=True #set true to run the speciesNAmeLists lists through the mass plot production process False runs a single plot
nplot = 7 #list to plot


switch=1
ice = True
papersize = "A6"
xaslog='linear'
withfactorlegendP1 = True
legendoutP1 = False
withfactorlegendS = True
legendoutS = False

imgparams=imgsize[papersize]

columnpath = "../VaryFromSolar/outputfiles"+str(switch)+"/"
if ice:
    plotspath = "../VaryFromSolar/"+papersize+xaslog+"facSepPlots"+str(switch)
else:
    plotspath = "../VaryFromSolar/"+papersize+xaslog+"facNoIcePlots"+str(switch)

if os.path.isdir(plotspath) is not True:
    os.mkdir(plotspath)

pf.plt.rcParams['xtick.labelsize']=imgparams[5]
pf.plt.rcParams['ytick.labelsize']=imgparams[5]
#fig,axes=pf.plt.subplots(len(elements),len(varyfactor),figsize=(16,9))
#fig,axes=pf.plt.subplots(figsize=(16,9))#len(elements),len(varyfactor),figsize=(16,9))
#axes=axes.flatten()
#i=0

#for m , speciesNames in enumerate(speciesNameLists):
def plotchem(speciesNames):
    p=0
    speclist = ""
    factorlines = [mlines.Line2D([],[],color='w',label="factor",linestyle='None')]
    for j , strfactor in enumerate(varyfactorstr):
        factorlines.append( mlines.Line2D([],[],color='k',label=strfactor,linestyle=Linestles[j]))
    for l, s in enumerate(speciesNames):
        speclist += s + " "
    for k, e in enumerate(elements):
        i=0
        iprev = i
        specfactor = []
        title = "Varying " + e[0]
        #aiming to have 3 panes stacked vertically
        fig,axes=pf.plt.subplots(figsize=imgparams[0], num=p,clear=True)
        if legendoutP1:
            pf.plt.subplots_adjust(right=0.85,left=None,bottom=None,top=None)
        p=p+1
        #figcombo,axescombo=pf.plt.subplots(figsize=(16,9), sharey=True,num=p,clear=True)
        #p=p+1
#Separate out the phase 2 graph
        #figp2,axesp2=pf.plt.subplots(figsize=(16,9), sharey=True,num=p,clear=True)
        phase2startt=0
        phase2endt=0
        rtists = []
        for j , factor in enumerate(varyfactor):
            #pick species, any number is fine
            print("plotting: "+ speclist + e[0] +" " + str(factor))
            #title = ""
            #phase 1 data extract
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
          
            #phase 2 data notbeing used at mo so leave out
        #     if factor == 1 : 
        #         collfile = columnpath + "phase2-fullCtrl.dat"
        #         #lb=True
        #         #title = "Ctrl"
        #     else: 
        #         collfile = columnpath + "phase2-full"+e[0]+ "-" + str(factor).replace('.','_')+".dat"
        #         #lb=False
        #         #title = "Varying " + e[0]#+str(factor)
      		# #call read_uclchem. 
        #     time2,dens2,temp,abundances2=pf.read_uclchem(collfile,speciesNames)
        #   		#plot species returns the axis so we can further edit
        #     #stitch together the phase1 and phase2 series and plot
        #     timenp = np.asarray(time)
        #     timenp2 = np.asarray(time2)
        #     for l, s in enumerate(speciesNames):
        #         abundances[l].extend(abundances2[l])
            
        #     tnpp2 = timenp2 + timenp[-1]
        #     t = np.append(timenp,tnpp2)
            axis,rtist0=pf.plot_species(speciesNames,time,abundances,axes,ls=Linestles[j],lab=lb,lw=imgparams[4],ncol=6)#ax=axes[i])
            #axiscombo,rtist0=pf.plot_species(speciesNames,t,abundances,axescombo,ls=Linestles[j],lab=lb)#ax=axes[i])
            rtists.append(rtist0)#questionable use as plotspecies doesn't return the list of handles just the last ouneat the moment 
            
            # plot phase2 on its own
            #phase2startt=tnpp2[0]
            #phase2endt=tnpp2[-1]
      		#plot species and save to test.png, alternatively send dens instead of time.
            #axisp2,rtist2=pf.plot_species(speciesNames,tnpp2,abundances2,axesp2,ls=Linestles[j],lab=lb)#ax=axes[i])
        
        #here insert the factor linestyles at the end of the legned under a heading of factor
        if withfactorlegendP1:         
            lbls = ["species"]
            hdls = [mlines.Line2D([],[],color='w',label="species",linestyle='None')]
            hndls,labls = axes.get_legend_handles_labels()
            lbls.extend(labls)
            lbls.append("factor")
            lbls.extend(varyfactorstr)
            hdls.extend(hndls)
            hdls.extend(factorlines)
            
            
        axes.set(xscale=xaslog,yscale='log',ylim=(1e-18,1e-3),xlim=(1e0,time[-1]))#timenp[-1]
        axes.set_xlabel('Time / Years',fontsize=imgparams[2])
        if xaslog == 'linear':
            axes.ticklabel_format(axis='x',useMathText=True)
        axes.set_ylabel('X/H',fontsize=imgparams[2])
        axes.set_title(title+" phase1",fontsize=imgparams[3])
        
        if withfactorlegendP1 and legendoutP1:
            axes.legend(hdls,lbls,loc='upper left',bbox_to_anchor=(1.0,1.0),fontsize=imgparams[5],handlelength=1.5,borderaxespad=0.0)#loc='upper left',bbox_to_anchor=(0.0,0.0,1.0,1.0),
        elif withfactorlegendP1:
            axes.legend(hdls,lbls,loc='best',fontsize=imgparams[1])
        else:
            axes.legend(loc='best',fontsize=imgparams[1])
        
        fig.savefig(plotspath+"/phase1plot"+e[0]+"_"+speciesNames[0]+".png",dpi=300)
        
        # axes.set(xscale=xaslog,yscale='log',ylim=(1e-18,1e-3),xlim=(1e0,6.5e6))#t[-1]))
        # if xaslog == 'linear':
        #     axes.ticklabel_format(axis='x',useMathText=True)
        # axes.set_title(title+" phase1&2",fontsize=imgparams[3])
        # fig.savefig(plotspath+"/comboplot"+e[0]+"_"+speciesNames[0]+".png",dpi=300)

        # axes.set(xscale=xaslog,yscale='log',ylim=(1e-18,1e-3),xlim=(phase2startt,phase2startt+0.3e6))
        # if xaslog == 'linear':
        #     axes.ticklabel_format(axis='x',useMathText=True)
        # axes.set_title(title+" phase2",fontsize=imgparams[3])
        # fig.savefig(plotspath+"/phase2plot"+e[0]+"_"+speciesNames[0]+".png",dpi=300)

        

    p=p+1
    for k, e in enumerate(elements):
        i=0
        iprev = i
        specfactor = []
        #aiming to have 3 panes stacked vertically
        fig,axes=pf.plt.subplots(figsize=imgparams[0],num=p,clear=True)
        if legendoutS:
            pf.plt.subplots_adjust(right=0.85,left=None,bottom=None,top=None)
        for j , factor in enumerate(varyfactor):
            print("plotting static: "+speclist+e[0]+" "+str(factor))
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
            axis,rtist=pf.plot_species(speciesNames,time3,abundances3,axes,ls=Linestles[j],lab=lb,lw=imgparams[4],ncol=6)#ax=axes[i])
          
          		#plot species returns the axis so we can further edit
            if xaslog == 'linear':
                axis.set(xscale=xaslog,yscale='log',ylim=(1e-18,1e-3),xlim=(1e0,time3[-1]))
                axis.ticklabel_format(axis='x',useMathText=True)
            else:
                axis.set(xscale=xaslog,yscale='log',ylim=(1e-18,1e-3),xlim=(1e0,time3[-1]))
            axis.set_ylabel('X/H',fontsize=imgparams[2])
            axis.set_xlabel('Time / Years',fontsize=imgparams[2])
            
            
            axis.set_title(title+" static cloud",fontsize=imgparams[3])
            #axis.set_title(title)
            #here insert the factor linestyles at the end of the legned under a heading of factor
            if withfactorlegendS:         
                lbls = ["species"]
                hdls = [mlines.Line2D([],[],color='w',label="species",linestyle='None')]
                hndls,labls = axes.get_legend_handles_labels()
                lbls.extend(labls)
                lbls.append("factor")
                lbls.extend(varyfactorstr)
                hdls.extend(hndls)
                hdls.extend(factorlines)
            
            if withfactorlegendS and legendoutS:
                axes.legend(hdls,lbls,loc='upper left',bbox_to_anchor=(1.0,1.0),fontsize=imgparams[5],handlelength=1.5,borderaxespad=0.0)#loc='upper left',bbox_to_anchor=(0.0,0.0,1.0,1.0),
            elif withfactorlegendS:
                axes.legend(hdls,lbls,loc='best',fontsize=imgparams[1])
            else:
                axes.legend(loc='best',fontsize=imgparams[1])
            i=i+1
    
        #axes[0].text(.02,0.98,e[0],horizontalalignment="left",verticalalignment="top",transform=axes[0].transAxes)
        #the single plot per page version
#        axes.text(.02,0.98,e[0],horizontalalignment="left",verticalalignment="top",transform=axes.transAxes)
        
        #axes[3].text(.02,0.98,"Your Row",horizontalalignment="left",verticalalignment="top",transform=axes[3].transAxes)
        fig.savefig(plotspath+"/staticplot"+e[0]+"_"+speciesNames[0]+".png",dpi=300)
        pf.plt.clf()


if bulk:
    pool = mp.Pool(12)
    if ice:
        pool.map(plotchem, speciesNameLists)
    else:
        pool.map(plotchem, speciesNoiceNameLists)
    pool.close()
    pool.join()
else:
    if ice:
        plotchem(speciesNameLists[nplot])
    else:
        plotchem(speciesNoiceNameLists[nplot])

