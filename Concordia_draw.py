#!/usr/bin/env python
# encoding: utf-8
"""
Concordia_draw.py

Created by Mauricio Ibanez-Mejia on 2013-12-26.
Copyright (c) 2013 __MyCompanyName__. All rights reserved.
"""

import math
import numpy
#import matplotlib
import pylab as P
import matplotlib.pyplot as plt
from matplotlib.transforms import offset_copy
from matplotlib import cm
import sys
from mpl_toolkits.mplot3d import Axes3D

import CalculateConcordia as CC




#ratio75, error75, ratio68, error68, rho = numpy.loadtxt("/Users/jcibanezm/codes/Practice/Python/data_abs_two.txt", delimiter='\t', unpack=True)#.tolist()
#ratio75, error75, ratio68, error68, rho = numpy.loadtxt("/Users/Ibanez/Documents/Python_scripts/data_abs.txt", delimiter='\t', unpack=True)#.tolist()

#ratio75, error75, ratio68, error68, rho = numpy.loadtxt("data_abs.txt", delimiter='\t', unpack=True)#.tolist()
#ratio75, error75, ratio68, error68, rho = numpy.loadtxt("Fioretti.txt", delimiter='\t', unpack=True)#.tolist()
#ratio75, error75, ratio68, error68, rho = numpy.loadtxt("10MIGU21.txt", delimiter='\t', unpack=True)#.tolist()
#ratio75, error75, ratio68, error68, rho = numpy.loadtxt("n1000_4.txt", delimiter='\t', unpack=True)#.tolist()
#ratio75, error75, ratio68, error68, rho = numpy.loadtxt("CP40conc_4.txt", delimiter='\t', unpack=True)#.tolist()
#ratio75, error75, ratio68, error68, rho = numpy.loadtxt("CP40conc_4_PDF.txt", delimiter='\t', unpack=True)#.tolist()
#ratio75, error75, ratio68, error68, rho = numpy.loadtxt("cp40conc_4_less2200.txt", delimiter='\t', unpack=True)#.tolist()
#ratio75, error75, ratio68, error68, rho = numpy.loadtxt("CP40conc_11.txt", delimiter='\t', unpack=True)#.tolist()

# print 'Calculating min and max values'
# Determine the limits of the region to be studied.
max_values75 = ratio75 + (sigma+0.01) * error75
min_values75 = ratio75 - (sigma+0.01) * error75
max_limit75 = numpy.amax(max_values75)
min_limit75 = numpy.amin(min_values75)

max_values68 = ratio68 + (sigma+0.05) * error68
min_values68 = ratio68 - (sigma+0.05) * error68
max_limit68 = numpy.amax(max_values68)
min_limit68 = numpy.amin(min_values68)


#agestart = 10.
agestart = CC.age68(min_limit68)		#numpy.amin(ages68) - numpy.amin(ages68)*0.35
ageend   = CC.age75(max_limit75)		#numpy.amax(ages75) + numpy.amax(ages75)*0.35

ratio68start = min_values68		#val68fromage(agestart)
ratio68end   = max_limit68	#val68fromage(ageend)
ratio75start = min_limit75		#val75fromage(agestart)
ratio75end   = max_limit75		#val75fromage(ageend)

step = (ageend - agestart)/450.

if agestart < 100:
	frac_labels = 10.
	pow_labels = 1
elif agestart >= 100:
	frac_labels = 100.
	pow_labels = 1

if agestart-ageend > 1000:
	frac_labels=1000.
	pow_labels = 10

PlotLabelsRange = int((ageend-agestart)/frac_labels)*pow_labels+100
AgeLabelStep = 100
AgeLabelStep_logaxes = 100

print ''
print '*******************************************************************'
print "		Start the Concordia Calculations."
print '*******************************************************************'
print ''

################################    Start CALCULATE CONCORDIA    ###################################################

# Why 450 ???
# ------------- Start - Create arrays and calculate concordia ratios
agearray    = numpy.zeros(450)
conc68array = numpy.zeros(450)
conc75array = numpy.zeros(450)
conc67array = numpy.zeros(450)

#for i in range(450):
#	agearray.append(0)
#	conc68array.append(0)
#	conc75array.append(0)
#	conc67array.append(0)

for i in range(450):	# ------------- Calculates age data array
	agearray[i]=agestart+i*step+1e-10

for i in range(450):	# ------------- Calculates isotope ratio data array for concordia
	conc68array[i]=math.exp(lambda238*agearray[i])-1
	conc75array[i]=math.exp(lambda235*agearray[i])-1
	conc67array[i]=(U85)*((math.exp(lambda238*agearray[i])-1)/(math.exp(lambda235*agearray[i])-1))

# ------------- End - Create arrays and calculate concordia ratios


# Generate arrays to plot the Age labels in the concordia diagram
agelabelarray 	= []
labels68		= []
labels75		= []
labels67		= []
agelabelarray_logaxes 	= []
labels68logaxes			= []
labels75logaxes			= []
labels67logaxes			= []

for i in range(PlotLabelsRange):
	agelabelarray.append(0)
	labels68.append(0)
	labels75.append(0)
	labels67.append(0)
	agelabelarray_logaxes.append(0)
	labels68logaxes.append(0)
	labels75logaxes.append(0)
	labels67logaxes.append(0)


agelabelarray[0]=1.0e-10
agelabelarray_logaxes[0]=1.0e-10

for i in range(PlotLabelsRange-1):	# ------------- Calculates age data array
	agelabelarray[i+1]=int(agestart/100.)*100+(i)*AgeLabelStep+1e-10
	agelabelarray_logaxes[i+1] = int(agestart/100.)*100+(i)*AgeLabelStep_logaxes+1e-10

for i in range(PlotLabelsRange):	# ------------- Calculates isotope ratio data array for concordia
	labels68[i]=math.exp(lambda238*agelabelarray[i])-1
	labels75[i]=math.exp(lambda235*agelabelarray[i])-1
	labels67[i]=(U85)*((math.exp(lambda238*agelabelarray[i])-1)/(math.exp(lambda235*agelabelarray[i])-1))
	labels68logaxes[i]=math.exp(lambda238*agelabelarray_logaxes[i])-1
	labels75logaxes[i]=math.exp(lambda235*agelabelarray_logaxes[i])-1
	labels67logaxes[i]=(U85)*((math.exp(lambda238*agelabelarray_logaxes[i])-1)/(math.exp(lambda235*agelabelarray_logaxes[i])-1))

################################ End CALCULATE CONCORDIA ###################################################





######################################### 		Perform calculations 		#########################################

number_spots = len(ratio75) #Count number of analyses

max75ratio = numpy.amax(ratio75)	#Maximum and minimum values for 6/8 and 7/5 and their uncertainties in the imported array
min75ratio = numpy.amin(ratio75)
max75error = numpy.amax(error75)
min75error = numpy.amin(error75)
max68ratio = numpy.amax(ratio68)
min68ratio = numpy.amin(ratio68)
max68error = numpy.amax(error68)
min68error = numpy.amin(error68)

minxstep = 2*sigma*min75error/steps_x_axis	#Calculate the minimum step size necessary for the entire grid based on the most precise analysis
minystep = 2*sigma*min68error/steps_y_axis

xsteps_finalgrid = int((max_limit75 - min_limit75)/minxstep) #int((ratio75end - ratio75start)/minxstep)	#Calculate total number of steps for the final grid
ysteps_finalgrid = int((max_limit68 - min_limit68)/minystep)

print ''
print '---------------------------------------------------------------------'
print 'Properties for the Global grid.'
print ''
print '	Age Start	', agestart, ' Myr'
print '	Age End		', ageend, ' Myr'
print ''
print '	Min 207/235	', min_limit75, ' Max 207/235	', max_limit75
print '	Step-size	', minxstep, '  # steps	', xsteps_finalgrid
print ''
print '	Min 206/238	', min_limit68, ' Max 206/238	', max_limit68
print '	Step-size	', minystep, '  # steps	', ysteps_finalgrid
print ''
print 'Total # of cells =	', xsteps_finalgrid*ysteps_finalgrid
print '---------------------------------------------------------------------'
print ''


#print 'Generate general array'
x_finalgrid = numpy.linspace(min_limit75, max_limit75, xsteps_finalgrid, endpoint=True)	#Create arrays for X and Y for the final grid
y_finalgrid = numpy.linspace(max_limit68, min_limit68, ysteps_finalgrid, endpoint=True)
X, Y 		= numpy.meshgrid(x_finalgrid, y_finalgrid)	#Create final grid (individual analyses will be added to this one)

def initialize(X, Y):
	# why ?
	a = X*0 + Y*0 - 0.0001
	return a

#print 'Initialize general array'
Final_PDF = initialize(X, Y)

if old_concordia == 'yes':
	Temporal_PDF = Final_PDF

print ''
print ' ===================================================================='

xsteps_mem = xsteps_finalgrid*sys.getsizeof(x_finalgrid[0])/1.0e6
ysteps_mem = ysteps_finalgrid*sys.getsizeof(y_finalgrid[0])/1.0e6
cells = xsteps_finalgrid*ysteps_finalgrid
X_mem = cells*sys.getsizeof(X[0][0])/1.0e6
Y_mem = cells*sys.getsizeof(Y[0][0])/1.0e6
PDF_mem = cells*sys.getsizeof(Final_PDF[0][0])/1.0e6

print 'Sizes of arrays in MB'
print ''
print '	Size of 1D x final grid ', xsteps_mem
print '	Size of 1D y final grid ', ysteps_mem
print '	Size of 2D X, Y Grids 	', X_mem,',', Y_mem
print '	Size of the Final PDF	', PDF_mem
print '	Total Memory used (MB) =', xsteps_mem+ysteps_mem+2*X_mem+PDF_mem
print ' ===================================================================='
print ''

replot_concordia = 'yes'
if old_concordia == 'yes':
# Plot the traditional concordia diagram.
	fig = P.figure(figsize=(12,8.5))
	ax = P.subplot(1,1,1)

print '		Start calculating the individual PDFs'

while replot_concordia == 'yes':
	for i in range(number_spots):

		# Create a 0's array to find the position of the individual grid in the Ginal PDF
		x_node_diff = numpy.zeros(xsteps_finalgrid)
		y_node_diff = numpy.zeros(ysteps_finalgrid)

		# Measure the distance between the center of the PDF and the nodes in the general Grid.
		# Locate the closest node with the minimum distance.
		x_node_diff = abs(x_finalgrid - ratio75[i]).tolist()
		x_close_node_val = numpy.amin(x_node_diff).tolist()
		index_x_close_node = x_node_diff.index(x_close_node_val)
		y_node_diff = abs(y_finalgrid - ratio68[i]).tolist()
		y_close_node_val = numpy.amin(y_node_diff).tolist()
		index_y_close_node = y_node_diff.index(y_close_node_val)

		# Calculate the number of grid points of the individual PDF.
		pointsxgrid = int((2*sigma*error75[i])/minxstep)
		pointsygrid = int((2*sigma*error68[i])/minystep)

		# Calculate the half size of the PDF
		edgexgrid = int((sigma*error75[i])/minxstep)
		edgeygrid = int((sigma*error68[i])/minystep)

		# Locate the corners of the individual PDF to be matched with the Final PDF.
		xSW = x_finalgrid[index_x_close_node - edgexgrid]
		ySW = y_finalgrid[index_y_close_node + edgeygrid]
		xNE = x_finalgrid[index_x_close_node + edgexgrid - 1]
		yNE = y_finalgrid[index_y_close_node - edgeygrid - 1]

		# Calculate the control indexes in the Final PDF to match the individual PDF.
		controlindex_x_finalgrid = index_x_close_node - pointsxgrid/2	#Calculate a control node on the final_grid that should match the NW corner of the current spot grid
		controlindex_y_finalgrid = index_y_close_node - pointsygrid/2
		control_x_finalgrid = x_finalgrid[controlindex_x_finalgrid]
		control_y_finalgrid = y_finalgrid[controlindex_y_finalgrid]

		# End points of the individual PDF box.
		end_x_value = xSW + pointsxgrid * minxstep
		end_y_value = yNE - pointsygrid * minystep
		xdist = end_x_value - xSW
		ydist = yNE - end_y_value

		# Generate the corresponding array of the individial PDF
		x = numpy.linspace(xSW, end_x_value, pointsxgrid)#.tolist()
		y = numpy.linspace(yNE, end_y_value, pointsygrid)#.tolist()
		xx, yy = numpy.meshgrid(x,y)

		###############################################################################
		#						This is the important calculations
		###############################################################################
		# Calculate the correlated bivariate probability density function using eq 1 of Sircombe 06 paper.
		# And Normalize the PDF to 100.
		PDF = numpy.meshgrid(x,y, indexing='ij')
		# JCIM: I don't remember what is rho!!
		PDF = F(xx, yy, ratio75[i], ratio68[i], error75[i], error68[i], rho[i])
		# mmmmm. why do we normalize olf concordia?
		normalization_old_conc 	= numpy.max(PDF)
		normalization_value 	= numpy.sum(PDF)
		# JCIM: Ok I'm a bit confused here, why do we normalize stuff two times ???
		normPDF 				= PDF/(normalization_value)
		normPDF_old_conc = PDF/(normalization_old_conc)
		###############################################################################
		#
		#	print 'sum norm PDF', numpy.sum(normPDF)

		figure_range_x = pointsxgrid - 1
		figure_range_y = pointsygrid - 1

		# Print in Screen some properties of the individual PDF.
		print ''
		print '------------------------------------------------------------------------------------------'
		print 'Properties of the individual PDF'
		print ''
		print '	Sample number ', i
		print '	length of 2D array (x,y)	',len(x),',', len(y)
		print '	Normalization Value		', normalization_value
		print '	Max Probability			', F(ratio75[i], ratio68[i], ratio75[i], ratio68[i], error75[i], error68[i], rho[i])#/normalization_value
		print ''
		print '	SW corner index of individual '
		print '	PDF in the final grid (x,y)	',controlindex_x_finalgrid,',', controlindex_y_finalgrid
		print ''
		print '	SW corner value of the '
		print '	 individual PDF  (7/5,6/8)	',xSW,',', ySW
		print ''
		print '	NE corner value of the '
		print '	 individual PDF (7/5,6/8)	',xNE,',', yNE
		print ''
		print '-------------------------------------------------------------------------------------------'
		print ''

		print '	Fill the global PDF'

	    	for k in range(figure_range_x):
	   			for l in range(figure_range_y):
   					Final_PDF[controlindex_y_finalgrid + l][controlindex_x_finalgrid + k] = Final_PDF[controlindex_y_finalgrid + l][controlindex_x_finalgrid + k] + normPDF[l][k]
   					if old_concordia == 'yes':
						Temporal_PDF[controlindex_y_finalgrid + l][controlindex_x_finalgrid + k] = normPDF_old_conc[l][k]

		if old_concordia == 'yes':
			# Normal Color Contours
			line_contours = [0.05]
	#		line_contours = [0.32]
			plt.contour(X, Y, Temporal_PDF, line_contours, colors='red', linewidth=1)
			Temporal_PDF[:][:] = 0

	if old_concordia == 'yes':
		transOffset = offset_copy(ax.transData, fig=fig, x=-0.05, y=0.08, units='inches')
		plt.plot(conc75array, conc68array)

		xdown 	= min_limit75
		xup		= max_limit75
		ydown	= min_limit68
		yup		= max_limit68
		age_min = labels_ratio_limit(xdown, labels75)
		age_max = labels_ratio_limit(xup  , labels75)

		logaxes = raw_input("Do you want logarithmic axes ?")

		jkp = 0
		for i in range(age_min, age_max):
			if logaxes == 'no':
				plt.plot((labels75[i],),(labels68[i],), 'k+')
				plt.text(labels75[i], labels68[i], '%d' % agelabelarray[i], transform=transOffset)
			else:
				if agelabelarray[i] < 1000:
					plt.plot((labels75[i],),(labels68[i],), 'k+')
					plt.text(labels75[i], labels68[i], '%d' % agelabelarray[i], transform=transOffset)
				elif agelabelarray[i] > 1000:
					plt.plot((labels75[i+jkp],),(labels68[i+jkp],), 'k+')
					plt.text(labels75[i+jkp], labels68[i+jkp], '%d' % agelabelarray[i+jkp], transform=transOffset)
					jkp=jkp+1

		ax.set_xlim(xdown, xup)
		ax.set_ylim(ydown, yup)

		if logaxes == 'yes':
			ax.set_xscale('log')
			ax.set_yscale('log')
		plt.show()

		if logaxes == 'yes':
			fig.savefig("Conc_95_log.pdf", format='PDF')
		else:
			fig.savefig("Conc_95_norm.pdf", format='PDF')


	if old_concordia == 'yes':
		replot_concordia = raw_input("Do you want to replot the old concordia ?")
	else:
		replot_concordia = 'no'


# Normalize the Final PDF.
final_norm = numpy.max(Final_PDF)
Final_PDF = Final_PDF/final_norm*100
print ''



################################ PLOTTING ###################################################

#output = raw_input("What type of plot do you want, '2D' or '3D': ")
#print output

#print logaxes

replot = 'yes'

### Limits of the box in the plot.
#xdown 	= min_limit75
#xup		= max_limit75
#ydown	= min_limit68
#yup		= max_limit68

xdown 	= 6.0e-2
xup		= max_limit75
ydown	= 1.0e-2
yup		= max_limit68

#xup		= proper_limit(7	, x_finalgrid)
#yup		= proper_limit(0.4	, y_finalgrid)

age_min = labels_ratio_limit(xdown, labels75)
age_max = labels_ratio_limit(xup  , labels75)

print '	Plot the Concordia diagram with the cummulative PDFs'

while replot == 'yes':

	logaxes = raw_input("Do you want logarithmic axes ? ('yes' or 'no') ")
	show_plot = raw_input("Do you want the plot to be shown ? ")

	if output == '2D':
		######### 2D Plot ##################
		fig = P.figure(figsize=(12,8.5))
		ax = P.subplot(1,1,1)

		transOffset = offset_copy(ax.transData, fig=fig, x=-0.25, y=0.08, units='inches')

#		if logplot == 'no':
		# Normal Color Contours
		color_contours = [1,10,20,30,40,50,60,70,80,90,99,100]
		line_contours = [1,25,50,75,100]
		line2_contours = [1,10,20,30,40,50,60,70,80,90,100]

#		conc = plt.contourf(X,Y, Final_PDF, color_contours, cmap='jet')
		conc = plt.contourf(X,Y, Final_PDF, color_contours, cmap='jet')
		plt.contour(X, Y, Final_PDF, line2_contours, colors='black', linewidth=.2)

		##############    TO DO: allow for logarithmic contours.    #################
#		if logplot == 'yes':
#			# Logarithmic Color contours
#			log_color_contours 	= [1,10,20,30,40,50,60,70,80,90,99,100]
#			log_line_contours 	= [1,25,50,75,100]
#			log_line2_contours 	= [1,10,20,30,40,50,60,70,80,90,100]
#
#			conc = plt.contourf(X,Y, numpy.log(Final_PDF,10), color_contours, cmap='hot_r')
#			plt.contour(X, Y, numpy.log(Final_PDF,10), line2_contours, colors='black', linewidth=.5)

		plt.plot(conc75array, conc68array)

#		for i in range(PlotLabelsRange):
		jkp = 0
		jkp2= 0
		for i in range(age_min, age_max):
			#plt.plot((labels75[i],),(labels68[i],), 'ro')
#			plt.plot((labels75[i],),(labels68[i],), 'k+')
			if logaxes == 'no':
				plt.plot((labels75[i],),(labels68[i],), 'k+')
				plt.text(labels75[i], labels68[i], '%d' % agelabelarray[i], transform=transOffset)
			else:
				if agelabelarray[i] < 1000:
					plt.text(labels75[i], labels68[i], '%d' % agelabelarray[i], transform=transOffset)
					plt.plot((labels75[i],),(labels68[i],), 'k+')
				elif agelabelarray[i] > 1000:
					plt.plot((labels75[i+jkp],),(labels68[i+jkp],), 'k+')
					plt.text(labels75[i+jkp], labels68[i+jkp], '%d' % agelabelarray[i+jkp], transform=transOffset)
					jkp=jkp+1

		fig.colorbar(conc, shrink=0.8, aspect=10)

		ax.set_xlim(xdown, xup)
		ax.set_ylim(ydown, yup)

		if logaxes == 'yes':
			ax.set_xscale('log')
			ax.set_yscale('log')

		if show_plot == 'yes':
			plt.show()


	if output == '3D':
		color_contours = [1,10,20,30,40,50,60,70,80,90,99,100]
		############# 3D Plot ################
		fig = P.figure(figsize=(8,5.5))
		ax = P.subplot(1,1,1, projection='3d')

		transOffset = offset_copy(ax.transData, fig=fig, x=-0.05, y=0.08, units='inches')

		surface = ax.plot_surface(X,Y,Final_PDF, rstride=1, cstride=1, cmap='hot_r', linewidth=0, antialiased=False)
		fig.colorbar(surface, shrink=1.0, aspect=5)

		ax.set_xlabel('207/235')
		ax.set_ylabel('206/238')

	#	ax.set_xlim3d(min75ratio, max75ratio);
	#	ax.set_ylim3d(min68ratio, min68ratio);
		ax.set_zlim3d(1,100);

		plt.plot(conc75array, conc68array, 'r')
		for i in range(PlotLabelsRange):
			plt.plot((labels75[i],),(labels68[i],),(0,), 'ro')
	#		plt.text(labels75[i], labels68[i], '%d' % agelabelarray[i], transform=transOffset)
	## Have not been able to put the age labels.


	#	surface = ax.plot_surface(X,Y,Final_PDF, rstride=1, cstride=1, cmap='hot_r', linewidth=0, antialiased=False)
		if show_plot == 'yes':
			plt.show()


	# Ask the User for the
	save_plot = 'no'
	save_plot = raw_input("Do you want to save the plot (yes or no): ")
	print save_plot

	if save_plot == 'yes':
		filenamehere = raw_input("Enter the name of the plot: ")
		fig.savefig(filenamehere+".pdf", format='PDF')

	replot = raw_input("Do you want do re-do the plot (yes or no)? ")

	if replot == 'yes':
		print ''
		update_box_limits = raw_input("Do you want to update the box limits ")
		if update_box_limits == 'yes':
			print "Remember x-limits have to be between ", xdown, "and", xup
			print "Remember y-limits have to be between ", ydown, "and", yup
			xdown 	= input("low limit in x = ")
			xup		= input("up  limit in x = ")
			ydown	= input("low limit in y = ")
			yup		= input("up  limit in y = ")
			xdown 	= proper_limit(xdown, x_finalgrid)
			xup		= proper_limit(xup	, x_finalgrid)
			ydown	= proper_limit(ydown, y_finalgrid)
			yup		= proper_limit(yup	, y_finalgrid)

		age_min = labels_ratio_limit(xdown, labels75)
		age_max = labels_ratio_limit(xup  , labels75)


################################ Values to Print ###################################################
#print conc68array
#print ratio75[2]
#print number_spots
#
#print max75ratio
#print min75ratio
#print max75error
#print min75error
#print max68ratio
#print min68ratio
#print max68error
#print min68error

#print minxstep
#print minystep

#print agestart
#print ageend
#
#print ratio68start
#print ratio68end
#print ratio75start
#print ratio75end

#print xsteps_finalgrid
#print ysteps_finalgrid

#print X, Y

#print x_node_diff
#print x_close_node_val
#print index_x_close_node
#print index_y_close_node

#print xSW
#print xNE
#print ySW
#print yNE
#
#print controlindex_x_finalgrid
#print controlindex_y_finalgrid
#print control_x_finalgrid
#print control_y_finalgrid
