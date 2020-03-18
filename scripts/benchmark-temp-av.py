#! /usr/bin/python
#
import codecs
import math
import os
import time
import numpy as np
import matplotlib
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import FormatStrFormatter, LogLocator, MultipleLocator

# Specify the path for the ASCII data files
path = '/home/jon/Documents/zupcx4/UCL_CHEMS/uclpdr/'

# Specify the prefix and suffix for the model output data files
prefix = 'test'
suffix = '.out'


# No user changes needed beyond this point...
# Convert a string representation of a number to a float, handling the
# case where the exponent is larger than +/-99 and the 'E' is missing
def number(string):
    try:
        return float(string)
    except ValueError:
        if string.count('+') > 0: fixedString = string[:string[2:].index('+')+2] + 'E' + string[string[2:].index('+')+2:]
        if string.count('-') > 0: fixedString = string[:string[2:].index('-')+2] + 'E' + string[string[2:].index('-')+2:]
        return float(fixedString)

# Read the physical cloud properties
def read_cloud_properties(inputfile):
    inputfile.seek(0)
    data = inputfile.readline().split()
    gasDensity = [] ; dustDensity = [] ; gasTemperature = [] ; dustTemperature = [] ; fuvFlux = [] ; xrayFlux = []
    data = inputfile.readline().split()
    while data != []:
        gasDensity.append(number(data[1]))
        dustDensity.append(number(data[2]))
        gasTemperature.append(number(data[3]))
        dustTemperature.append(number(data[4]))
        fuvFlux.append(number(data[5]))
        xrayFlux.append(number(data[6]))
        data = inputfile.readline().split()
    nParticles = len(gasDensity)
    return nParticles, gasDensity, dustDensity, gasTemperature, dustTemperature, fuvFlux, xrayFlux

# Read the particle coordinates and the visual extinction along each ray to the PDR surface
def read_visual_extinctions(inputfile):
    inputfile.seek(0)
    data = inputfile.readline().split()
    particleCoordinates = []; visualExtinction = []
    data = inputfile.readline().split()
    while data != []:
        nRays = len(data)-4
        particleCoordinates.append([number(data[i]) for i in range(1,4)])
        visualExtinction.append([number(data[i]) for i in range(4,len(data))])
        data = inputfile.readline().split()
    return nRays, particleCoordinates, visualExtinction



# Read the various data files
print('\nReading cloud properties...')
dataFilename = prefix+'.prop.out'
input = codecs.open(path+dataFilename, encoding='utf-8', mode='r')
nParticles, gasDensity, dustDensity, gasTemperature, dustTemperature, fuvFlux, xrayFlux = read_cloud_properties(input)
input.close()

print('Reading visual extinctions...')
dataFilename = prefix+'.av.out'
input = codecs.open(path+dataFilename, encoding='utf-8', mode='r')
nRays, particleCoordinates, visualExtinction = read_visual_extinctions(input)
input.close()


print(len(gasTemperature),len(visualExtinction))
fig,ax=plt.subplots(figsize=(16,9))
ax.plot(visualExtinction,gasTemperature,label="UCLPDR")
ax.set(xscale='log',yscale='log')


model_no,cloud_size,av=np.loadtxt("output/av_grid.dat",unpack=True,delimiter=",",skiprows=1)

avs=[]
temps=[]
for model in model_no:
    data=np.loadtxt(f"output/av_grid/{model:.0f}.dat")
    avs.append(np.mean(data[-20:-1,3]))
    temps.append(np.mean(data[-20:-1,2]))

ax.plot(avs,temps,label="UCLCHEM")
ax.legend()
fig.savefig("output/benchmark.png")