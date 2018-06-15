#!/usr/bin/env python
# encoding: utf-8
"""
Concordia_draw.py

Created by Mauricio Ibanez-Mejia & Juan Ibanez-Mejia  2013-2018.
Copyright (c) 2013 __Ibanez-Mejia__. All rights reserved.
"""

print ("==================================================================")
print ("Running Concordia_draw.py")
print ("==================================================================")

import math
import numpy as np
#import matplotlib
import pylab as P
import matplotlib.pyplot as plt
from matplotlib.transforms import offset_copy
from matplotlib.colors import LogNorm
from mpl_toolkits.mplot3d import Axes3D
import sys
from timeit import time

import CalculateConcordia as CC
import PlotEllipse as PE

startTime = time.time()
# Decaying constants.
# JCIM: Do we need this here?
lambda235 = 9.8485E-4
lambda238 = 1.55125E-4
U85       = 137.818

#desired level of confidence and minimum number of steps for the most precise analysis
sigma = 4

#"Read the input parameters"
# Define the type of plot, 2D PDF of the concordia or a 3D surface plot, and set the axis to normal or logarithm.
output        = '2D'
log_axes      = False
filenamehere  = 'test'
old_concordia = False
PDFres        = 20
numConcSteps  = 100
#"Read the data"

data_dir   = "./"
out_dir    = "./"
input_data = ""

f = open("concordia_input_parameters.txt", 'r')

# Run this while not end of file ...
#for i in range(11):
print("")
print("Input parameters:")
for entry in f:
#    entry = f.readline()
    if entry.split() != []:
        if entry.split()[0] != '#':
            par_name  = entry.split()[0]
            par_input = entry.split()[2]
            #print(par_name, par_input)
            exec("%s = %s" %(par_name, par_input))
            print("\t%s \t\t= %s" %(par_name, par_input))

f.close()

dictD = {"info":"Dictionary of the data."}
dictD["in"], dictD["out"], dictD["res"] = "%s/%s"%(data_dir, input_data), "%s"%out_dir, PDFres

# Check for the input parameters:
# Make a Function to check for the necessary input parameters...
CC.TestInputParameters(dictD)

ratio75, error75, ratio68, error68, rho = np.loadtxt("%s/%s"%(data_dir,input_data), delimiter='\t', unpack=True)#.tolist()

number_spots = len(ratio75) #Count number of analyses

agearray, conc68array, conc75array, conc67array = CC.StartConcordiaLine(ratio68, error68, ratio75, error75, sigma)
agelabelarray, labels68, labels75, labels67     = CC.ComputeConcordiaLabels(agearray, old_concordia=True)

CC.TestInputParameters(dictD, memory=True, ratio68=ratio68, error68=error68, ratio75=ratio75, error75=error75, sigma=sigma)

#print("Labels age",agelabelarray)
#print("Labels 68", labels68)
#print("Labels 75", labels75)

#print("Labels 67", labels67)

###################################################################
#                  Plot the Old Concordia
###################################################################

if old_concordia == True:

    #Start plotting the figure
    fig = plt.figure(figsize=(12, 8.5))
    ax = fig.add_subplot(111)
    ax.plot(conc75array, conc68array, "-k")
    plot_kwargs = {'color':'r','linestyle':'-','linewidth':2,'alpha':0.8}
    for i in range(number_spots):
        # I need to revisit the ellipse plotting algorithm because of the difference range of the axes.
        PE.plot_ellipse(semimaj=error75[i]/100., semimin=error68[i],phi=rho[i], ax=ax, x_cent=ratio75[i], y_cent=ratio68[i], plot_kwargs=plot_kwargs)

    # ToDo: Deal with the labels being plotted outside the plot.
    for i in range(len(agelabelarray)):
        ax.plot(labels75[i],labels68[i], 'k+')
        ax.text(labels75[i], labels68[i], '%d' % agelabelarray[i])

    ax.set_xlim(0.9*np.min(conc75array), 1.1*np.max(conc75array))
    ax.set_ylim(0.9*np.min(conc68array), 1.1*np.max(conc68array))

    ax.tick_params(axis='both', which='major', length=7, width=2,  labelsize=12, direction="in")
    ax.tick_params(axis='both', which='minor', length=4, width=1.5, labelsize=12, direction="in")

    #if log_axes == True:
    #    ax.set_xscale("log")
    #    ax.set_yscale("log")

    fig.savefig("%s/%s_standard.pdf"%(out_dir,figure_name), format="pdf")

else:
    ###################################################################
    #                  Plot the New density PDF Concordia
    ###################################################################
    # Add option for logarithmic axes
    # Add option for logarithmic density distribution.

    X, Y, PDF = CC.ComputePDF(ratio75, error75, ratio68, error68, rho, PDFres=PDFres, sigma=sigma, log_axes=log_axes)

    fig = plt.figure(figsize=(12, 8.5))
    ax = fig.add_subplot(111)

    ax.plot(conc75array, conc68array, "-k", linewidth=0.5)

    for i in range(len(agelabelarray)):
        if labels75[i] > np.min(X) and labels75[i] <= np.max(conc75array):
            ax.plot(labels75[i],labels68[i], 'k+', color="r", linewidth=0.5)
            ax.text(labels75[i]+0.01*labels75[i], labels68[i]-0.01*labels68[i], '%d' % agelabelarray[i])

    if log_contour == False:
        color_contours = [1,10,20,30,40,50,60,70,80,90,100]
        line_contours = [1,25,50,75,100]
        conc_bar = ax.contourf(X,Y,PDF, color_contours, cmap='RdYlBu_r', alpha=0.9)
        ax.contour(X, Y, PDF, line_contours, colors='black', linewidths=0.1)
        cbar = fig.colorbar(conc_bar)
    else:
        log_contours   = [-4,-3.5,-3,-2.5,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2]
        log_line_contours = [-4,-2,0,2]
        conc_bar = ax.contourf(X,Y,np.log10(PDF+1.0e-4), log_contour, cmap='RdYlBu_r', alpha=0.9)
        ax.contour(X, Y, np.log10(PDF), log_line_contours, colors='black', linewidths=0.1)
        cbar = fig.colorbar(conc_bar)

    ax.set_xlim(np.min(X), np.max(X))
    ax.set_ylim(np.min(Y), np.max(Y))

    ax.tick_params(axis='both', which='major', length=7, width=2,  labelsize=12, direction="in")
    ax.tick_params(axis='both', which='minor', length=4, width=1.5, labelsize=12, direction="in")

    if log_axes == True:
        ax.set_xscale("log")
        ax.set_yscale("log")

    fig.savefig("%s/%s_PDF.pdf"%(out_dir,figure_name), format="pdf")



endTime = time.time()

totalTime = endTime - startTime
print("")
print("---------------- Done ----------------------")

if totalTime < 60:
    print("Time taken in the calculation was %.2f sec"%totalTime)
elif totalTime > 60 and totalTime < 3600:
    print("Time taken in the calculation was %.2f min"%(totalTime/60.))
elif totalTime < 3600.:
    print("Time taken in the calculation was %.2f hours"%(totalTime/3600.))
#print("Time taken in calculation: %.2f"%(endTime-startTime))
print("")
