import numpy as np

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
    a = 1/(2*math.pi*xerr*yerr*math.sqrt(1-rho**2))
    b = -0.5*1/(1-rho**2)
    c = ((xx - xmu)/xerr)**2 - 2*rho*(((xx - xmu)*(yy - ymu))/(xerr*yerr))+((yy - ymu)/yerr)**2
    return a*np.exp(b*c)


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

def AddConcordiaLabels(agearray, old_concordia=False, axislog=False, ax=None, data_out=False):
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

    ax : matplotlib axis property
        A pre-created matplotlib axis

    data_out : boolean
        A flag to return the axis labels array.
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

    PlotLabelsRange = int((ageend-agestart)/frac_labels)*pow_labels+100
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


    for i in range(PlotLabelsRange):	# ------------- Calculates isotope ratio data array for concordia
    	labels68[i]=math.exp(lambda238*agelabelarray[i])-1
    	labels75[i]=math.exp(lambda235*agelabelarray[i])-1
    	labels67[i]=(U85)*((math.exp(lambda238*agelabelarray[i])-1)/(math.exp(lambda235*agelabelarray[i])-1))
    	labels68logaxes[i]=math.exp(lambda238*agelabelarray_logaxes[i])-1
    	labels75logaxes[i]=math.exp(lambda235*agelabelarray_logaxes[i])-1
    	labels67logaxes[i]=(U85)*((math.exp(lambda238*agelabelarray_logaxes[i])-1)/(math.exp(lambda235*agelabelarray_logaxes[i])-1))

    ax.blablabla
