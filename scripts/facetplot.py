#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 20 21:10:09 2020

@author: norman
"""
#Plot in an table form to splurge a set of similar graphs on a page
#it reads full UCLCHEM output and saves a plot of the abudances of select species 
import os
import multiprocessing as mp
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import plotfunctions as pf
#Element ids and initial values
# This part can be substituted with any choice of grid
elements=[["C","fc",2.6e-04], ["N","fn",6.1e-05], ["O","fo",4.6e-04], ["S","fs",1.318e-05],
           ["Mg","fmg",3.981e-05], ["Si","fsi",1.0e-07],["Cl","fcl",3.162e-07]]#,["F","ff",3.6e-08]]

imgsize = {"A7":[(3.5,2),3.0,4.0,4.0,0.5,4],
           "A6":[(5.2,3.5),'xx-small','x-small','small',1.0,7],
           "A5":[(7.4,3.2),'xx-small','x-small','small',1.0,6],
           "A4":[(10,6.8),'small','medium','large',2.0,9]}

#list the species we want graphing from the current data output of the model
#speciesNames = ("SO SO2 S2 N NH NH2 NH3 HNC NO NO2 OCN HNCO HCS O2 H2O").split()
speciesNameLists = []
speciesNameLists.append(["S1",("SO SO2 S2 #SO2 #SO").split(),("SO SO2 S2 SO2ice SOice").split()])
speciesNameLists.append(["S2",("HS H2S #HS #H2S").split(),("HS H2S HSice H2Sice").split()])
speciesNameLists.append(["S3",("HCS OCS H2S2 #OCS #H2S2").split(),("HCS OCS H2S2 OCSice H2S2ice").split()])
speciesNameLists.append(["N1",("NH NH2 #NH #NH2").split(),("NH NH2 NHice NH2ice").split()])
speciesNameLists.append(["N2",("NH3 HNC #NH3 #HNC").split(),("NH3 HNC NH3ice HNCice").split()])
speciesNameLists.append(["N3",("NO NO2 #NO #NO2").split(),("NO NO2 NOice NO2ice").split()])
speciesNameLists.append(["N4",("N HNCO #N #HNCO").split(),("N HNCO Nice HNCOice").split()])
speciesNameLists.append(["N5",("HCN H2CN #HCN #H2CN").split(),("HCN H2CN HCNice H2CNice").split()])
speciesNameLists.append(["O1",("H2O HNO #H2O #HNO").split(),("H2O HNO H2Oice HNOice").split()])
speciesNameLists.append(["O2",("O O2 OH #O #O2 #OH").split(),("O O2 OH Oice O2ice OHice").split()])
speciesNameLists.append(["C1",("C CH CH2 #C #CH #CH2").split(),("C CH CH2 Cice CHice CH2ice").split()])# #C #CH #CH2
speciesNameLists.append(["C2",("CH3 CH4 #CH3 #CH4").split(),("CH3 CH4 CH3ice CH4ice").split()])
speciesNameLists.append(["C3",("C3H2 CH3CCH #C3H2").split(),("C3H2 CH3CCH C3H2ice").split()])
speciesNameLists.append(["C4",("CO CO2 #CO #CO2").split(),("CO CO2 COice CO2ice").split()])
speciesNameLists.append(["Si1",("SIC SIH2 #SIC #SIH4").split(),("SIC SIH2 SICice SIH4ice").split()])
speciesNameLists.append(["Mg",("CL HCL MG #HCL").split(),("Cl HCl Mg HClice").split()])
speciesNameLists.append(["Si2",("SIS SIC4 H2SIO SIO").split(),("SiS SiC4 H2SiO SiO").split()])
speciesNameLists.append(["e1",("MG MG+ E-").split(),("Mg Mg+ e-").split()])
speciesNameLists.append(["e2",("C C+ E-").split(),("C C+ e-").split()])
speciesNameLists.append(["e3",("S S+ E-").split(),("S S+ e-").split()])
speciesNameLists.append(["e4",("O O+ E- O2").split(),("O O+ e- O2").split()])
speciesNameLists.append(["e5",("C C+ E- MG MG+ S S+").split(),("C C+ e- Mg Mg+ S S+").split()])
speciesNameLists.append(["e6",("C C+ E- MG MG+ O O2 S S+").split(),("C C+ e- Mg Mg+ O O2 S S+").split()])

speciesNoiceNameLists = []
speciesNoiceNameLists.append(["S1",("SO SO2 S2").split(),("SO SO2 S2").split()])
speciesNoiceNameLists.append(["S2",("HS H2S").split(),("HS H2S").split()])
speciesNoiceNameLists.append(["S3",("HCS OCS H2S2").split(),("HCS OCS H2S2").split()])
speciesNoiceNameLists.append(["N1",("NH NH2").split(),("NH NH2").split()])
speciesNoiceNameLists.append(["N2",("NH3 HNC").split(),("NH3 HNC").split()])
speciesNoiceNameLists.append(["N3",("NO NO2").split(),("NO NO2").split()])
speciesNoiceNameLists.append(["N4",("N HNCO").split(),("N HNCO").split()])
speciesNoiceNameLists.append(["N5",("HCN H2CN").split(),("HCN H2CN").split()])
speciesNoiceNameLists.append(["O1",("H2O HNO").split(),("H2O HNO").split()])
speciesNoiceNameLists.append(["O2",("O O2 OH").split(),("O O2 OH").split()])
speciesNoiceNameLists.append(["C1",("C CH CH2").split(),("C CH CH2").split()])# #C #CH #CH2
speciesNoiceNameLists.append(["C2",("CH3 CH4").split(),("CH3 CH4").split()])
speciesNoiceNameLists.append(["C3",("C3H2 CH3CCH").split(),("C3H2 CH3CCH").split()])
speciesNoiceNameLists.append(["C4",("CO CO2").split(),("CO CO2").split()])
speciesNoiceNameLists.append(["Si1",("SIC SIH2").split(),("SIC SIH2").split()])
speciesNoiceNameLists.append(["Mg",("CL HCL MG ").split(),("Cl HCl Mg").split()])
speciesNoiceNameLists.append(["Si2",("SIS SIC4 H2SIO SIO").split(),("SiS SiC4 H2SiO SiO").split()])
speciesNoiceNameLists.append(["e1",("MG MG+ E-").split(),("Mg Mg+ e-").split()])
speciesNoiceNameLists.append(["e2",("C C+ E-").split(),("C C+ e-").split()])
speciesNoiceNameLists.append(["e3",("S S+ E-").split(),("S S+ e-").split()])
speciesNoiceNameLists.append(["e4",("O O+ E- O2").split(),("O O+ e- O2").split()])
speciesNoiceNameLists.append(["e5",("C C+ E- MG MG+ S S+").split(),("C C+ e- Mg Mg+ S S+").split()])
speciesNoiceNameLists.append(["e6",("C C+ E- MG MG+ O O2 S S+").split(),("C C+ e- Mg Mg+ O O2 S S+").split()])

varyfactor = [0.25, 0.5, 1, 2, 4]
#Linestles = [(0,(3,10,1,10)),(0,(3,5,1,5,1,5)),(0,()),(0,(3,5,1,5)),(0,(3,1,1,1))]
#Linestles = [(0,(4,2,1,2,1,2,1,2,1,2)),(0,(4,2,1,2,1,2)),(0,()),(0,(1,2,4,2,4,2)),(0,(1,2,4,2,4,2,4,2,4,2))]
Linestles = [(0,(1,2)),(0,(1,1)),(0,()),(0,(2,1,1,1)),(0,(2,1))]
Linestyles = [(1,2),(1,1),(),(2,1,1,1),(2,1)]

bulk=True #set true to run the speciesNoiceNameLists lists through the mass plot production process False runs a single plot
nplot = 17 #list to plot

modelname="UCLCHEM"

switch=1
ice = True
papersize = "A5"
xaslog='linear'
sns.set()

sns.set_context("paper")
#sns.axes_style(xscale=xaslog,yscale='log')
sns.set_palette("bright")
sns.set_style('whitegrid')

imgparams=imgsize[papersize]

columnpath = "../VaryFromSolar"+modelname+"/outputfiles"+str(switch)+"/"
if ice:
    plotspath = "../VaryFromSolar"+modelname+"/"+papersize+xaslog+"FGSepPlots"+str(switch)
else:
    plotspath = "../VaryFromSolar"+modelname+"/"+papersize+xaslog+"FGNoIcePlots"+str(switch)

if os.path.isdir(plotspath) is not True:
    os.mkdir(plotspath)

plt.rcParams['xtick.labelsize']=imgparams[5]
plt.rcParams['ytick.labelsize']=imgparams[5]

#plt.subplots(figsize=imgparams[0])
#fig,axes=pf.plt.subplots(len(elements),len(varyfactor),figsize=(16,9))
#fig,axes=pf.plt.subplots(figsize=(16,9))#len(elements),len(varyfactor),figsize=(16,9))
#axes=axes.flatten()
#i=0

#for m , speciesNames in enumerate(speciesNameLists):
def plotchem(speciesGroup):
    #speciesGroup=speciesNameLists[0]#use this to just do one element for a test
    speciesNames = speciesGroup[1]
    speciesDisplay = speciesGroup[2]
    groupname = speciesGroup[0]
    p=0
    for k, e in enumerate(elements):
        i=0
        iprev = i
        specfactor = []
        time = []
        abundances = []
        species = []
        abundscale = []
        varying=[]#willl be a different plot "varying" may be able to put in a grid of plots
        model = []
        timelim = 5.4e6
        title = modelname + " Varying " + e[0]
        #aiming to have 3 panes stacked vertically
        fig,axes=plt.subplots(figsize=imgparams[0], num=p,clear=True)
        plt.rcParams["axes.labelsize"]=imgparams[5]
        plt.rcParams["axes.titlesize"]="large"
        plt.rcParams["axes.formatter.use_mathtext"]=True
        plt.rcParams["font.size"]=8.0
        plt.rcParams["legend.fontsize"]="small"
        plt.rcParams["legend.borderaxespad"]=0.0
        plt.rcParams["legend.handletextpad"]=0.1
        #plt.rcParams["figure.autolayout"]=True

        p=p+1
        #figcombo,axescombo=pf.plt.subplots(figsize=(16,9), sharey=True,num=p,clear=True)
        #p=p+1
    #Separate out the phase 2 graph
        #figp2,axesp2=pf.plt.subplots(figsize=(16,9), sharey=True,num=p,clear=True)
        
        
        timeadded = False
    #        phase2startt=0
    #        phase2endt=0
        #iterate over all the multipliers held in varyfactor
        for j , factor in enumerate(varyfactor):
            #pick species, any number is fine
            #title = ""
            lb=False
            
            #phase 1 plot
            if factor == 1 : 
                collfile = columnpath + "phase1-fullCtrl.dat"
                #title = "Ctrl"
                lb=True
            else: 
                collfile = columnpath + "phase1-full"+e[0]+ "-" + str(factor).replace('.','_')+".dat"
                #title = "Varying " + e[0]#+str(factor)
                lb=False
      		#call read_uclchem. 
            t,dens,temp,abund=pf.read_uclchem(collfile,speciesNames)
            if t[-1] > timelim:
                timelim = t[-1]
            for l, s in enumerate(speciesDisplay):
                if len(abund[l]) != len(t):
                    print("Species "+s+"no values a=" +str(len(abund))+" t="+str(len(t)))
                    continue
                abundances.extend(abund[l])
                time.extend(t)
                species.extend([s]*len(t)) #.replace('ice','*')
                abundscale.extend([str(factor)]*len(t))
                varying.extend([e[0]]*len(t))#varying element means we could construct a page of plots
                model.extend(["collapse"]*len(t))
            
            #static model plot
            if factor == 1 : 
                collfile = columnpath + "static-fullCtrl.dat"
                title = "Ctrl"
                lb = True
            else: 
                collfile = columnpath + "static-full"+e[0]+ "-" + str(factor).replace('.','_')+".dat"
                title = modelname + "Varying " + e[0]#+str(factor)
                lb = False
      		#call read_uclchem. 
            t,dens,temp,abund=pf.read_uclchem(collfile,speciesNames)
            
            for l, s in enumerate(speciesDisplay):
                if len(abund[l]) != len(t):
                    print("Species "+s+"no values a=" +str(len(abund))+" t="+str(len(t)))
                    continue
                abundances.extend(abund[l])
                time.extend(t)
                species.extend([s]*len(t)) #.replace('ice','*')
                abundscale.extend([str(factor)]*len(t))
                varying.extend([e[0]]*len(t))#varying element means we could construct a page of plots
                model.extend(["static"]*len(t))
            # timenp = np.asarray(time)
            # if not timeadded:
            #     p1df["time"] = time
            #     timeadded = True
            # for l, s in enumerate(speciesNames):
            #     colname=s+str(factor)
            #     colname = s.replace('#','ice')
                
            #     specfactor.append(colname)
                
            #     p1df = pd.concat([p1df,pd.DataFrame({colname:abundances[l]})],axis=1)#column for each species with 
            
            #
      		#plot species and save to test.png, alternatively send dens instead of time.
            #axis,rtist0=pf.plot_species(speciesNames,time,abundances,axes,ls=Linestles[j],lab=lb)#ax=axes[i])
            #axis,rtist1=pf.plot_species(speciesNames,time,abundances,axes[1],ls=Linestles[j],lab=False)#ax=axes[i])
        #p1 = p1df.transpose()
        p1df = pd.DataFrame({"time":time,"abundances":abundances,"species":species,"factor":abundscale,"varying":varying,"model":model})
        
        g = sns.FacetGrid(p1df,col="model",row="varying",aspect = 1.0)#,legend_out=True)
        
        
        g.map_dataframe(sns.lineplot, x="time",y="abundances",hue="species",style="factor",dashes=Linestyles,linewidth=imgparams[4],ci=None) #data=p1df,legend="brief",ax=axes,
        g.set(xscale=xaslog,yscale='log',ylim=(1e-18,5e-3),xlim=(1e0,timelim)) 
        g.set_axis_labels('Time years','X/H')
        g.add_legend(fontsize=imgparams[5])#,loc='upper left',bbox_to_anchor=(1.0,1.0))
        g.fig.suptitle("Abundances "+modelname, fontsize=8)
        g.fig.subplots_adjust(top=0.9,left=0.1,right=0.9)

        #axes.set(xscale=xaslog,yscale='log',ylim=(1e-18,1e-3),xlim=(1e0,t[-1]))
        #axes.set_xlabel('Time / Years',fontsize=imgparams[2])
        #if xaslog == 'linear':
        #    axes.ticklabel_format(axis='x',useMathText=True)
       # axes.set_ylabel('X/H',fontsize=imgparams[2])
        #axes.set_title(title+" static cloud",fontsize=imgparams[3])
        #axes.legend(loc='best',fontsize=imgparams[1])
        plt.savefig(plotspath+"/facetplot"+e[0]+"_"+groupname+".png",dpi=300)  
             
            
    
        #axes[0].text(.02,0.98,e[0],horizontalalignment="left",verticalalignment="top",transform=axes[0].transAxes)
        #the single plot per page version
#        axes.text(.02,0.98,e[0],horizontalalignment="left",verticalalignment="top",transform=axes.transAxes)
        
        #axes[3].text(.02,0.98,"Your Row",horizontalalignment="left",verticalalignment="top",transform=axes[3].transAxes)
#        fig.savefig(plotspath+"/staticplot"+e[0]+"_"+speciesNames[0]+".png",dpi=300)
        #pf.plt.clf()


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

