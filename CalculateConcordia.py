import numpy as np
import math
import sys

# Decaying constants.
lambda235 = 9.8485E-4
lambda238 = 1.55125E-4
U85       = 137.818

# Functins used to find the age of a sample given the Pb206/U238  & Pb207/U235 Ratio
def age68(ratio68):
    """
    Document ...
    """
    return np.log(ratio68+1)/lambda238

def age75(ratio75):
    """
    Document ...
    """
    return np.log(ratio75+1)/lambda235

# Functins used to find the ratios of Pb206/U238  & Pb207/U235 given the age.
def val68fromage(age68):
    """
    Document ...
    """
    return np.exp(lambda238*age68)-1

def val75fromage(age75):
    """
    Document ...
    """
    return np.exp(lambda235*age75)-1

# What do this functions do ??
def proper_limit(limit_in, array):
    """
    Document ...
    """
    aux = 1
    for i in range(len(array)):
        diff = abs(limit_in - array[i])
    if diff < aux:
        aux = diff
        limit = array[i]
    return limit

def labels_ratio_limit(ratio_limit_in, array):
    """
    Document ...
    """
    aux = 10
    index = len(array)
    for i in range(len(array)):
        diff = abs(ratio_limit_in - array[i])
    if diff < aux:
        aux = diff
        index = i
    return index

# Equation 1 in Sircombe 06 Paper. correlated bivariate normal distribution.
def F(xx, yy, xmu, ymu, xerr, yerr, rho):	#Define function for the 2-d bivariate correlated PDF
    """
    Document ...
    """
    a = 1./(2.*math.pi*xerr*yerr*np.sqrt(1.-rho**2.))
    b = -0.5*1./(1.-rho**2.)
    c = ((xx - xmu)/xerr)**2. - 2.*rho*(((xx - xmu)*(yy - ymu))/(xerr*yerr))+((yy - ymu)/yerr)**2
    return a*np.exp(b*c)


def TestInputParameters(dictD, memory=False, ratio68=None, error68=None, ratio75=None, error75=None, sigma=None):
    """
    This function tests the input data and output directory to see if they exist.
    If they don't. print a Warning on screen and kill the code.
    """
    if memory == False:
        try:
           fn=open(dictD["in"],"r")
        except IOError:
           print("")
           print("\tError: the file %s does not appear to exist."%dictD["in"])
           print("\tCheck that you have the correct 'data_dir' and 'input_data' in the 'concordia_input_parameters.txt' file")
           print("")
           sys.exit()
    else:
        max_values75 = ratio75 + (sigma+0.01) * error75
        min_values75 = ratio75 - (sigma+0.01) * error75
        max_limit75  = np.amax(max_values75)
        min_limit75  = np.amin(min_values75)
        max_values68 = ratio68 + (sigma+0.05) * error68
        min_values68 = ratio68 - (sigma+0.05) * error68
        max_limit68  = np.amax(max_values68)
        min_limit68  = np.amin(min_values68)
        minxstep = 2*sigma*np.min(error75)/dictD["res"]	#Calculate the minimum step size necessary for the entire grid based on the most precise analysis
        minystep = 2*sigma*np.min(error68)/dictD["res"]
        xsteps_finalgrid = int((max_limit75 - min_limit75)/minxstep) #int((ratio75end - ratio75start)/minxstep)	#Calculate total number of steps for the final grid
        ysteps_finalgrid = int((max_limit68 - min_limit68)/minystep)

        entry = float(1.0)
        cells = xsteps_finalgrid*ysteps_finalgrid
        X_mem = cells*sys.getsizeof(entry)/1.0e6

        total_memory = 3*X_mem

        print ("")
        print (' ====================================================================')
        print ('                        Memory usage of the calculation:')
        print ('')
        print ('	Total Memory required (MB) = %.1f'%(total_memory))
        print (' ====================================================================')
        print ('')

        # Ask the system what is my available RAM. If I'm using more than half, ask if
        # I really want to run with this resolution.

    return 0

def StartConcordiaLine(ratio68, error68, ratio75, error75, sigma, numConcSteps=450):
    """
    Document ...
    """
    ages68 = age68(ratio68)
    ages75 = age75(ratio75)

    max_values75 = ratio75 + (sigma+0.01) * error75
    min_values75 = ratio75 - (sigma+0.01) * error75
    max_limit75  = np.amax(max_values75)
    min_limit75  = np.amin(min_values75)

    max_values68 = ratio68 + (sigma+0.05) * error68
    min_values68 = ratio68 - (sigma+0.05) * error68
    max_limit68  = np.amax(max_values68)
    min_limit68  = np.amin(min_values68)

    agestart = age68(min_limit68)		#numpy.amin(ages68) - numpy.amin(ages68)*0.35
    ageend   = age75(max_limit75)		#numpy.amax(ages75) + numpy.amax(ages75)*0.35

    ratio68start = min_values68		#val68fromage(agestart)
    ratio68end   = max_limit68	#val68fromage(ageend)
    ratio75start = min_limit75		#val75fromage(agestart)
    ratio75end   = max_limit75		#val75fromage(ageend)

    #"Start a figure and add the concordia line"
    AgeStep = (ageend - agestart)/numConcSteps

    if agestart < 100.:
    	frac_labels = 10.
    	pow_labels = 1
    elif agestart >= 100:
    	frac_labels = 100.
    	pow_labels = 1

    if agestart-ageend > 1000:
    	frac_labels=1000.
    	pow_labels = 10

    # I need to think about this
    PlotLabelsRange = int((ageend-agestart)/frac_labels)*pow_labels+100
    AgeLabelStep = 100
    AgeLabelStep_logaxes = 100
    #"Ask if you want to print new or old condordia"

    ################################    Start CALCULATE CONCORDIA    ###################################################
    # ------------- Start - Create arrays and calculate concordia ratios
    agearray    = np.zeros(numConcSteps)
    conc68array = np.zeros(numConcSteps)
    conc75array = np.zeros(numConcSteps)
    conc67array = np.zeros(numConcSteps)

    #print("Agestart = %.2f, numConcStep = %i, AgeStep = %.2g"%(agestart, numConcSteps, AgeStep))
    agearray = agestart + np.arange(numConcSteps)*AgeStep + 1e-10

    conc68array = np.exp(lambda238*agearray)-1.
    conc75array = np.exp(lambda235*agearray)-1.
    conc67array = U85*conc68array/conc75array

    return agearray, conc68array, conc75array, conc67array


def ComputeConcordiaLabels(agearray, old_concordia=False, axislog=False):
    """
    Add the labels to the concordia line depending on the range covered by the age array.
    Whether we plot the old or the new concordia in linear or logarithmic axis
    An easy to use function for plotting labels in the concordia line in Python 2.7!

    The function takes the age array and creates a label array depending on some properties of the plot.

    agearray : numpy array
        concordia age array with between the youngest and oldest point in the figure.

    old_concordia : boolean (True, False)
        Am I plotting the old concordia with the red ellipses or the new age PDF concordia.

    axislog : boolean
        what is the scale of the plot axis, linear or logarithmic.

    return:
        agelabelarray, labels68, labels67, labels75
    """

    agestart = agearray[0]
    ageend   = agearray[-1]

    if agestart < 100:
    	frac_labels = 10.
    	pow_labels = 1
    elif agestart >= 100:
    	frac_labels = 100.
    	pow_labels = 1

    if agestart-ageend > 1000:
    	frac_labels=1000.
    	pow_labels = 10

    PlotLabelsRange = int((ageend-agestart)/frac_labels)*pow_labels+10
    AgeLabelStep = 100
    AgeLabelStep_logaxes = 100

    # Generate arrays to plot the Age labels in the concordia diagram
    agelabelarray 	= np.zeros(PlotLabelsRange)
    labels68		= np.zeros(PlotLabelsRange)
    labels75		= np.zeros(PlotLabelsRange)
    labels67		= np.zeros(PlotLabelsRange)
    agelabelarray_logaxes 	= np.zeros(PlotLabelsRange)
    labels68logaxes			= np.zeros(PlotLabelsRange)
    labels75logaxes			= np.zeros(PlotLabelsRange)
    labels67logaxes			= np.zeros(PlotLabelsRange)

    agelabelarray[0]=1.0e-10
    agelabelarray_logaxes[0]=1.0e-10

    agelabelarray         = int(agestart/100.)*100+AgeLabelStep*np.arange(PlotLabelsRange)+1e-10
    agelabelarray_logaxes = int(agestart/100.)*100+AgeLabelStep_logaxes*np.arange(PlotLabelsRange)+1e-10

    labels68        = np.exp(lambda238*agelabelarray)-1.
    labels75        = np.exp(lambda235*agelabelarray)-1.
    labels67        = (U85)*((np.exp(lambda238*agelabelarray)-1.)/(np.exp(lambda235*agelabelarray)-1.))
    labels68logaxes = np.exp(lambda238*agelabelarray_logaxes)-1.
    labels75logaxes = np.exp(lambda235*agelabelarray_logaxes)-1.
    labels67logaxes = (U85)*((np.exp(lambda238*agelabelarray_logaxes)-1.)/(np.exp(lambda235*agelabelarray_logaxes)-1.))

    if axislog:
        agelabelarray = agelabelarray_logaxes
        labels68      = labels68logaxes
        labels67      = labels67logaxes
        labels75      = labels75logaxes

    return agelabelarray, labels68, labels75, labels67

def initialize(X, Y):
    a = X*0 + Y*0 - 0.0001
    return a

def ComputePDF(ratio75, error75, ratio68, error68, rho, PDFres=10, sigma=4, log_axes=False):
    """
    Add Documentation.
    """
    number_spots = len(ratio75) #Count number of analyses

    ages68   = age68(ratio68)
    ages75   = age75(ratio75)
    agestart = np.min(ages68)		#numpy.amin(ages68) - numpy.amin(ages68)*0.35
    ageend   = np.max(ages75)		#numpy.amax(ages75) + numpy.amax(ages75)*0.35

    max_values75 = ratio75 + (sigma+0.01) * error75
    min_values75 = ratio75 - (sigma+0.01) * error75
    max_limit75  = np.amax(max_values75)
    min_limit75  = np.amin(min_values75)

    max_values68 = ratio68 + (sigma+0.05) * error68
    min_values68 = ratio68 - (sigma+0.05) * error68
    max_limit68  = np.amax(max_values68)
    min_limit68  = np.amin(min_values68)

    max75ratio = np.amax(ratio75)	#Maximum and minimum values for 6/8 and 7/5 and their uncertainties in the imported array
    min75ratio = np.amin(ratio75)
    max75error = np.amax(error75)
    min75error = np.amin(error75)
    max68ratio = np.amax(ratio68)
    min68ratio = np.amin(ratio68)
    max68error = np.amax(error68)
    min68error = np.amin(error68)

    minxstep = 2*sigma*min75error/PDFres	#Calculate the minimum step size necessary for the entire grid based on the most precise analysis
    minystep = 2*sigma*min68error/PDFres

    xsteps_finalgrid = int((max_limit75 - min_limit75)/minxstep) #int((ratio75end - ratio75start)/minxstep)	#Calculate total number of steps for the final grid
    ysteps_finalgrid = int((max_limit68 - min_limit68)/minystep)

    print ('')
    print ('---------------------------------------------------------------------')
    print ('Properties of the Global grid.')
    print ('')
    print ('	Age Start	%.2f Myr'%agestart)
    print ('	Age End		%.2f Myr'%ageend)
    print ('')
    print ('	Min 207/235	%.2g Max 207/235  %.2g	'%(min_limit75, max_limit75))
    print ('	Step-size	%.2f # steps %i'%(minxstep, xsteps_finalgrid))
    print ('')
    print ('	Min 206/238   %.2g	 Max 206/238	%.2g'%(min_limit68, max_limit68))
    print ('	Step-size	%.2g   # steps	%i'%(minystep, ysteps_finalgrid))
    print ('')
    print ('    Total # of cells =	%i'%(xsteps_finalgrid*ysteps_finalgrid))
    print ('---------------------------------------------------------------------')
    print ('')

    # Do I want to have a different grid for logarithmic axis.
    #print 'Generate general array'
    if log_axes == False:
        x_finalgrid = np.linspace(min_limit75, max_limit75, xsteps_finalgrid, endpoint=True)	#Create arrays for X and Y for the final grid
        y_finalgrid = np.linspace(max_limit68, min_limit68, ysteps_finalgrid, endpoint=True)
    else:
        x_finalgrid = np.logspace(np.log10(min_limit75), np.log10(max_limit75), xsteps_finalgrid, endpoint=True)	#Create arrays for X and Y for the final grid
        y_finalgrid = np.logspace(np.log10(max_limit68), np.log10(min_limit68), ysteps_finalgrid, endpoint=True)

    # Create the basegrid for the calculations.
    X, Y 		= np.meshgrid(x_finalgrid, y_finalgrid)	#Create final grid (individual analyses will be added to this one)
    concordiaPDF = initialize(X, Y)

    for i in range(number_spots):
        concordiaPDF += F(X, Y, ratio75[i], ratio68[i], error75[i], error68[i], rho[i])

    # normalization
    normalization = np.max(concordiaPDF)
    #normalization = np.sum(concordiaPDF)
    concordiaPDF = concordiaPDF / normalization * 100.

    return X, Y, concordiaPDF

def CalculateMemoryUsage(X, Y):

    xsteps_finalgrid = np.shape(X)[1]
    ysteps_finalgrid = np.shape(X)[2]

    #xsteps_mem = xsteps_finalgrid*sys.getsizeof(x_finalgrid[0])/1.0e6
    #ysteps_mem = ysteps_finalgrid*sys.getsizeof(y_finalgrid[0])/1.0e6
    cells = xsteps_finalgrid*ysteps_finalgrid
    X_mem = cells*sys.getsizeof(X[0][0])/1.0e6
    Y_mem = cells*sys.getsizeof(Y[0][0])/1.0e6
    PDF_mem = cells*sys.getsizeof(X[0][0])/1.0e6

    print ('Sizes of arrays in MB')
    print ('')
    #print '	Size of 1D x final grid ', xsteps_mem
    #print '	Size of 1D y final grid ', ysteps_mem
    print ('	Size of 2D X = %.1g, Y =  %.1g grid'%(X_mem,Y_mem))
    print ('	Size of the Final PDF %.2g	'%PDF_mem)
    #print '	Total Memory used (MB) =', xsteps_mem+ysteps_mem+2*X_mem+PDF_mem
    print ('	Total Memory used (MB) = %.1f'%(2*X_mem+PDF_mem))
    print (' ====================================================================')
    print ('')
