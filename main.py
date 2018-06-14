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
PDFres        = 20
numConcSteps  = 450
#"Read the data"

print ('')
print ('*******************************************************************')
print ("		Start the Concordia Calculations.")
print ('*******************************************************************')
print ('')

ratio75, error75, ratio68, error68, rho = np.loadtxt("Sircombe.txt", delimiter='\t', unpack=True)#.tolist()

number_spots = len(ratio75) #Count number of analyses

agearray, conc68array, conc75array, conc67array = CC.StartConcordiaLine(ratio68, error68, ratio75, error75, sigma)

#Start plotting the figure
fig = plt.figure(figsize=(12, 8.5))
ax = fig.add_subplot(111)
ax.plot(conc75array, conc68array, "-k")
plot_kwargs = {'color':'r','linestyle':'-','linewidth':2,'alpha':0.8}
for i in range(number_spots):
    # I need to revisit the ellipse plotting algorithm because of the difference range of the axes.
    PE.plot_ellipse(semimaj=error75[i], semimin=error68[i],phi=rho[i], ax=ax, x_cent=ratio75[i], y_cent=ratio68[i], plot_kwargs=plot_kwargs)

CC.AddConcordiaLabels(agearray, old_concordia=True, ax=ax)

fig.savefig("test_oldConcordia.pdf", format="pdf")

#"Add concordia labels"
#    "I could make a figure that initialized this"


#if old_concordia == True:
#    "Run the Plot Ellipse code"

#else:
#    "Create the grid with the desired resolution"
#    "Loop over the analysis and fill the PDF"
#    "Plot the 2D concordia density distribution"
