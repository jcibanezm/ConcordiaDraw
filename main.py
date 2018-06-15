#!/usr/bin/env python
# encoding: utf-8
"""
Concordia_draw.py

Created by Mauricio Ibanez-Mejia on 2013-12-26.
Copyright (c) 2013 __MyCompanyName__. All rights reserved.
"""

import math
import numpy as np
#import matplotlib
import pylab as P
import matplotlib.pyplot as plt
from matplotlib.transforms import offset_copy
from matplotlib.colors import LogNorm
from mpl_toolkits.mplot3d import Axes3D
import sys

import CalculateConcordia as CC
import PlotEllipse as PE

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
logaxes       = 'yes'
filenamehere  = 'test'
old_concordia = 'no'
PDFres        = 60
numConcSteps  = 100
#"Read the data"

data_dir = "/Users/juan/Dropbox/codes/Geochronology/concordia_draw/Data"
out_dir  = "/Users/juan/Dropbox/codes/Geochronology/concordia_draw/Output"

print ('')
print ('*******************************************************************')
print ("		Start the Concordia Calculations.")
print ('*******************************************************************')
print ('')

ratio75, error75, ratio68, error68, rho = np.loadtxt("%s/%s"%(data_dir,"Sircombe.txt"), delimiter='\t', unpack=True)#.tolist()

number_spots = len(ratio75) #Count number of analyses

agearray, conc68array, conc75array, conc67array = CC.StartConcordiaLine(ratio68, error68, ratio75, error75, sigma)
agelabelarray, labels68, labels75, labels67     = CC.ComputeConcordiaLabels(agearray, old_concordia=True)

#print("Labels age",agelabelarray)
#print("Labels 68", labels68)
#print("Labels 75", labels75)

#print("Labels 67", labels67)

###################################################################
#                  Plot the Old Concordia
###################################################################

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

fig.savefig("%s/%s"%(out_dir,"test_oldConcordia.pdf"), format="pdf")

###################################################################
#                  Plot the New density PDF Concordia
###################################################################

X, Y, PDF = CC.ComputePDF(ratio75, error75, ratio68, error68, rho, PDFres=PDFres, sigma=sigma)

fig = plt.figure(figsize=(12, 8.5))
ax = fig.add_subplot(111)
ax.plot(conc75array, conc68array, "-k")

color_contours = [1,10,20,30,40,50,60,70,80,90,100]
line_contours = [1,25,50,75,100]
line2_contours = [1,10,20,30,40,50,60,70,80,90,100]

conc_bar = ax.contourf(X,Y,PDF, color_contours, cmap='RdYlBu_r', alpha=0.9)
ax.contour(X, Y, PDF, line2_contours, colors='black', linewidths=0.1)

#ax.set_xscale("log")
#ax.set_yscale("log")

cbar = fig.colorbar(conc_bar)
#cbar.ax.set_yticklabels(['< -1', '0', '> 1'])  # vertically oriented colorbar

#cax = ax.imshow(PDF, cmap="RdYlGn_r", aspect='auto')#,
#extent=[0.9*np.min(conc75array), 1.1*np.max(conc75array), 0.9*np.min(conc68array), 1.1*np.max(conc68array)])

# ToDo: Deal with the labels being plotted outside the plot.
for i in range(len(agelabelarray)):
    if labels75[i-1] <= np.max(conc75array):
        ax.plot(labels75[i],labels68[i], 'k+', color="r", linewidth=0.5)
        ax.text(labels75[i]+0.01*labels75[i], labels68[i]-0.01*labels68[i], '%d' % agelabelarray[i])

ax.set_xlim(0.95*np.min(conc75array), 1.05*np.max(conc75array))
ax.set_ylim(0.95*np.min(conc68array), 1.05*np.max(conc68array))

fig.savefig("%s/%s"%(out_dir,"test_newConcordia.pdf"), format="pdf")


#else:
#    "Create the grid with the desired resolution"
#    "Loop over the analysis and fill the PDF"
#    "Plot the 2D concordia density distribution"
