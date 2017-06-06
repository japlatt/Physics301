import numpy as np
import emcee
import kplr
from scipy.interpolate import interp1d
import batman



def Fw(t ,t0, F0, sigma = 0.01):
    """
    t - single time t
    Compute locally weighted F(t) given samples of t0 and

    f0, weighted by some timescale sigma

    returns 0 if no points within ~5 sigma.

    """
    weights = np.exp(-0.5*(t-t0)**2/sigma**2)/np.sqrt(2*np.pi*sigma**2)

    if np.isclose(np.sum(weights),0):
        F = 0.0
    else:
        F = np.average(F0,weights=weights)

    return F
    
def Query_database(exoplanet_name):
    """
    Query the Kepler database, and retrieve the lightcurve for a given
    exoplanet.  Clean the lightcurve and normalize to 1.  Return
    the time samples, flux samples, and the error on the flux samples.
    """
    client = kplr.API()
    planet = client.planet(exoplanet_name)
    lcs = planet.get_light_curves(short_cadence=False)
    time = np.zeros(0)
    flux = np.zeros(0)
    ferr = np.zeros(0)
    
    for lc in lcs:
        with lc.open() as f:
            hdu_data = f[1].data
            time = np.append(time,hdu_data["time"][np.isfinite(hdu_data["PDCSAP_FLUX"])])
            flux = np.append(flux,hdu_data["PDCSAP_FLUX"][np.isfinite(hdu_data["PDCSAP_FLUX"])]/np.nanmedian(hdu_data["PDCSAP_FLUX"]))
            ferr = np.append(ferr,hdu_data["PDCSAP_FLUX_ERR"][np.isfinite(hdu_data["PDCSAP_FLUX"])]/np.nanmedian(hdu_data["PDCSAP_FLUX"]))
    
    return time , flux , ferr


def Scan_for_transits(time,flux,lowerlim,upperlim,deltaT):
    """
    Scan a time series for the presence of transits occuring with
    times between lowerlim and upperlim, spaced by deltaT.  This is
    slow, so it should only be performed once.
    
    This function will return 2 arrays:
    
    depth:      The maximum signal depth.
    deltaT:     The time intervals used in the search.
    """
    # bins and flux allow us to enhance our signal by using multiple 
    # flux measurements made during different orbits to tease out the signal
    bins = np.linspace(0,1,300)
    Flux = np.zeros(len(bins)-1)
    
    deltaT = np.arange(lowerlim,upperlim,deltaT)
    depth = np.zeros(deltaT.shape)
    sigma = np.zeros(deltaT.shape)
    for j,t in enumerate(deltaT):
        print t
        for i in range(len(bins)-1):
            inds = np.logical_and(((time%t)/t>bins[i]),(time%t)/t<=bins[i+1])

            Flux[i] = np.median(flux[inds])
        depth[j] = np.min(Flux)
        sigma[j] = np.min((Flux-np.median(Flux))/np.std(Flux[np.where(abs(Flux-np.median(Flux))/np.std(Flux)<3)]))
        
    return depth ,deltaT , sigma
    
def Identify_transits(deltaT,depth,sigma):
    """
    parse the full array of deltaT and depth, and identify the individual transit candidates, including
    their time of peak signal and the depth of that signal.  Use the significance array to determine when a 
    signal has actually dipped.
    
    return cleaned list of transit candidates (any depth greater than 5 standard deviations)
    """
    
    i = 10
    
    transit_list = []
    
    while i<len(deltaT):
        if sigma[i] <= -5:
            # transit candidate, search for brightest nearby signal
            Transit_time = 0.0
            Transit_depth = np.inf
            significance = 0.0
            for j in range(i-10,i):
                if depth[j] < Transit_depth:
                    Transit_time = deltaT[j]
                    Transit_depth = depth[j]
                    significance = abs(sigma[j])
            # have scanned from behind the first detected location, now want to scan ahead
            # do so until the processed signal is back above 1 sigma
            j = 0
            while sigma[j+i] <= -5:
                if depth[j+i] < Transit_depth:
                    Transit_time = deltaT[j+i]
                    Transit_depth = depth[j+i]
                    significance = abs(sigma[j+i])
                j += 1
            transit_list.append([Transit_time,Transit_depth,significance])
            i += j
            
        else:
            i +=1
    
    Transit_array = np.zeros([len(transit_list),3])
    for i in range(len(transit_list)):
        Transit_array[i,0] = transit_list[i][0]
        Transit_array[i,1] = transit_list[i][1]
        Transit_array[i,2] = transit_list[i][2]
    return Transit_array
        
    
def Get_scan_info(deltaT,depth,sigma,Flux,Time):
    """
    Pass a list of transit candidates with period deltaT and transit depth depth.  For all of these
    calculate the following:
    
    Ntransits:          The number of transits that occur when the lightcurve is folded on this timescale
                        This is useful for discriminating real transits from harmonics.
    
    Transitdepthstd:    The standard deviation of the measured fluxes at the peak of the transit signal.  
                        We expect this should be larger if the timescale deltaT is a fraction of the true
                        transit signal, because it will place transits in resonance, but also non-transits
                        at the same location.
    
    Return an array containing the transit period, depth, Number of transits detected, and the sttdev of the 
    flux at peak depth.
    """
    depthsorted = [y for (y,x) in sorted(zip(depth,deltaT))]
    Timesorted  = [x for (y,x) in sorted(zip(depth,deltaT))]
    sigmasorted = [x for (y,x) in sorted(zip(depth,sigma))]
    
    data = np.zeros([len(deltaT),5])
    
    bins = np.linspace(0,1,300)
    flux = np.zeros(len(bins)-1)
    fluxstd = np.zeros(len(bins)-1)
    
    for j in range(len(deltaT)):
        t = Timesorted[j]
        print t
        Ncycles = 0
        for i in range(len(bins)-1):
            inds = np.logical_and(((Time%t)/t>bins[i]),(Time%t)/t<=bins[i+1])
            flux[i] = np.median(Flux[inds])
            fluxstd[i] = np.std(Flux[inds])
            
        for i in range(1,len(flux)):
            if (flux[i]-np.median(flux)<-4*np.std(flux)) & (flux[i-1]-np.median(flux)>=-4*np.std(flux)):
                Ncycles += 1
            
        data[j,0] = Timesorted[j]
        data[j,1] = depthsorted[j]
        data[j,2] = Ncycles
        data[j,3] = fluxstd[np.argmin(flux)] 
        data[j,4] = sigmasorted[j]
    
    return data
    
def Remove_harmonics(data):
    """
    Takes the output of the Get_scan_info function.  Removes any lines which have more than 1 transit, or which are
    a 1/int fraction of the strongest 5 transits (with only 1 detected dip).
    
    Returns the cleaned array
    """
    # First, clean away any candidates with multiple transits.
    first_cleaned_list = []
    for i in range(data.shape[0]):
        if np.isclose(data[i,2],1.0):
            first_cleaned_list.append(data[i])
    
    # second, remove any remaining candidates that are 1/int fractions of the strongest signals
    # iterate until the end of the list is reached.
    i=0
    while i < len(first_cleaned_list):
        for j in range(len(first_cleaned_list)-1,i):
            for k in range(2,10):
                if (abs(first_cleaned_list[j][0]/first_cleaned_list[i][0] - 1./float(k)) < 0.001):
                    first_cleaned_list.pop(j)
        i += 1
    cleaned_array = np.zeros([len(first_cleaned_list),len(first_cleaned_list[0])])
    for i in range(cleaned_array.shape[0]):
        cleaned_array[i,:] = first_cleaned_list[i]
    return cleaned_array 
        


class DEOptimizer(object):
    """
    This will be our program that we use to optimize our posterior, prior to
    sampling.  It utilizes Self_Adaptive Differential Evolution to explore
    the parameters while avoiding falling into local minima.
    """
    def __init__(self,NWalkers):
        """
        Initialize the optimizer
        """
        self.NWalkers = NWalkers
    
    def Initialization(self,priorup,priordn):
        """
        Set the initalization region.  Initialize Vectors
        with a uniform random distribution between priorup and
        priordn.
        """
        self.VECTORS = np.random.random([len(priorup),self.NWalkers])
        for i in range(len(priorup)):
                self.VECTORS[i,:] *= (priorup[i]-priordn[i])
                self.VECTORS[i,:] += priordn[i]
                
        self.DONORS = np.zeros(self.VECTORS.shape)
        self.TRIALS = np.zeros(self.VECTORS.shape)
        self.CURRENT_FX = np.ones(self.VECTORS.shape[1])*np.inf
        self.StrategyProbs = np.array([0.25,0.25,0.25,0.25])
        return
        
    def Mutation(self):
        """
        Mutate the vectors.
        """
        for i in range(self.VECTORS.shape[1]):
            options = range(self.VECTORS.shape[1])
            options.remove(i)
            choices = np.random.choice(options,3,replace=False)
            F = np.random.normal(0.5,0.3)
            self.DONORS[:,i] = self.VECTORS[:,choices[0]]+F*(self.VECTORS[:,choices[1]]-self.VECTORS[:,choices[2]])
        return
        
    def Recombination(self):
        """
        Recombine the donor vectors with the current vectors
        to create the trial vectors, including crossover.
        """
        for i in range(self.VECTORS.shape[1]):
            CR = np.random.normal(0.5,0.1)
            for j in range(self.VECTORS.shape[0]):
                if ((np.random.random() <= CR) or (j == np.random.choice(range(0,self.DONORS.shape[0])))):
                    self.TRIALS[j,i] = self.DONORS[j,i]
                else:
                    self.TRIALS[j,i] = self.VECTORS[j,i]
        return
        
    def Selection(self,ftominimize,x,y,yerr):
        """
        call ftominimize with the trial vectors, accept
        if it improves the fit, and replace existing vectors
        with the trial vectors.
        
        ftominimize should accept 4 arguments:  The x coords
        of datapoints, the y coords of sample points, the
        errors on the y coords, and an array of parameters,
        which the function handles internally to compute the model.
        """
        for i in range(self.VECTORS.shape[1]):
            f1 = self.CURRENT_FX[i]
            f2 = ftominimize(x,y,yerr,self.TRIALS[:,i])
            
            if f2 <= f1:
                self.VECTORS[:,i] = self.TRIALS[:,i]
                self.CURRENT_FX[i] = f2
            else:
                pass
        return
        
    def UpdateStrategyProbs(self):
        """
        Update the probabilities that a given strategy will be attempted, as well
        as CR and F for recombination.
        """
        return
    def Optimize(self,Niter,convThresh,function,x,y,yerr):
        """
        Run the optimization for either Niter iterations, or until the rms scatter 
        over all the dimensions is lower than convThresh
        """
        for i in range(Niter):
            self.Mutation()
            self.Recombination()
            self.Selection(function,x,y,yerr)
            
            if np.mean(np.std(self.VECTORS,axis=1)) < convThresh:
                break
                
        print "Optimization has Converged"
        print "Best Parameters:  " , self.VECTORS[:,np.argmin(self.CURRENT_FX)]
        return
def correlatefn(time, fluxfn, dt):
    """
    Compute auto-correlation function of signal

    @parameters:
    time - array of times

    flux - the signal to be autocorrelated over. Needs
    to be a continuous function from 0 - t, with masked out points
    set to 0.

    n_coors - Max number of data points to correlated over.  Will
    compute correlations from 0 - n_coors


    filt - Boolean to filter or not to filter

    std_filter - Mask out outliers with signal > std_filter*std

    """
    delt = time[1] - time[0]
    time = np.arange(time[0], time[-1], delt)

    coor = []
    flux = fluxfn(time)
    mean = np.mean(flux)
    standard_dev = np.std(flux)

    for delta in dt:
        time_mod = np.arange(time[0]+1.1*delta, time[-1]-1.1*delta, delt)
        I1 = fluxfn(time_mod)
        I2 = fluxfn(time_mod+delta)

        inds1 = (I1 >= 0)
        inds2 = (I2 >= 0)

        inds = inds1 & inds2
        N = sum(inds)
        coor.append(sum((I1[inds]*I2[inds]))/N)
    return np.array(coor)


def correlate(flux, n_coors = 80,filt = False, std_filter = 5):
    """
    Compute auto-correlation function of signal

    @parameters:
    flux - the signal to be autocorrelated over. Needs
    to be a continuous function from 0 - t, with masked out points
    set to 0.

    n_coors - Max number of data points to correlated over.  Will
    compute correlations from 0 - n_coors


    filt - Boolean to filter or not to filter

    std_filter - Mask out outliers with signal > std_filter*std

    """
    coor = np.zeros(n_coors)
    mean = np.nanmean(flux)
    standard_dev = np.nanstd(flux)
    for delta in xrange(n_coors):
        I1 = np.copy(flux)
        I2 = np.array(list(flux[delta:]) + list(flux[0:delta]))

        inds1 = (I1 >= 0.5)
        inds2 = (I2 >= 0.5)

        if(filt):
            inds1 = inds1 & (I1 >= mean - std_filter*standard_dev) & \
                            (I1 <=  mean + std_filter*standard_dev)
                
            inds2 = inds2 & (I2 >= mean - std_filter*standard_dev) & \
                            (I2 <=  mean + std_filter*standard_dev)

        inds = inds1 & inds2
        N = sum(inds)
        coor[delta] = sum((I1[inds]*I2[inds]))/N
    return coor
    
def Cn(time,flux,deltat):
    """
    Compute the correlation function for a given delta T scanned over all times
    
    time is list of time points.
    
    flux is scipy.interp1d for flux samples.
    
    deltaT is the time offset to scan through.
    """
    C = 0.0
    N = 0
    
    for i in range(len(time)):
        Fi = flux(time)
        Fipdt = Fw(time+deltat)
        if not np.isclose(Fipdt,0):
            C = np.sum( Fi*Fipdt )
            N = len(Fipdt) - np.sum(np.isclose(Fipdt,0))
    
    C /= float(N)
    return C
    
def loglikelihood(time,flux,fluxerr,parameters):
    """
    Compute the chi2 for a given set of parameters.
    """
    chi2 = np.sum((flux-TransitModel(time,parameters))**2/fluxerr**2)
    return chi2

def TransitModel(time,parameters):
    """
    Use the Batman transit modeling code to produce
    a transit model brightness (normalized to 1), at
    each time point.
    """
    
    params     = batman.TransitParams()
    params.t0  = parameters[0]                      # time of inferior conjunction
    params.per = parameters[1]                      # orbital period
    params.rp  = parameters[2]                      # planet radius (in units of stellar radii)
    params.a   = parameters[3]                      # semi-major axis (in units of stellar radii)
    params.inc = parameters[4]                      # orbital inclination (in degrees)
    params.ecc = parameters[5]                      # eccentricity
    params.w   = parameters[6]                      # longitude of periastron (in degrees)
    params.u   = [parameters[7],parameters[8]]      # limb darkening coefficients [u1, u2]
    params.limb_dark = "quadratic"                  # limb darkening model
    
    model = batman.TransitModel(params, time)
    return model.light_curve(params)

def Get_parameter_guesses(time,Flux,period_guess):
    """
    Fold the lightcurve on the period guess, and use the binning scheme to estimate the peak transit depth,
    transit width, and transit time.  These will give the planetary radius, orbital separation, and 
    transit time.  Return a vector with guesses of the parameters to feed to the DE optimizer (maybe
    upper and lower bounds as well).
    """
    time -= np.min(time)
    
    bins = np.linspace(0,period_guess,300)
    flux = np.zeros(len(bins)-1)
    
    for i in range(len(bins)-1):
        inds = np.logical_and(((time%period_guess)>bins[i]),(time%period_guess)<=bins[i+1])
        flux[i] = np.median(Flux[inds])
    t0 = bins[np.argmin(Flux)]+np.diff(bins)[0]/2.
    rp = np.sqrt(1.-np.min(flux))
    ap = (2*np.pi*(np.max(bins[np.where(flux < np.median(flux) - 0.8*(np.median(flux) -\
         np.min(flux)))]) - np.min(bins[np.where(flux < np.median(flux) - 0.8*(np.median(flux) - np.min(flux)))])))**-1. # Rough estimate
    inc = 90.
    ecc = 0.0
    w = 90.
    u1 = 0.3
    u2 = 0.3
    
    return np.array([t0,period_guess,rp,ap,inc,ecc,w,u1,u2])
    
def Refine_period_estimate(time,Flux,period_guess):
    """
    Zoom in further on the period estimate, and return it.
    """
    t,d,sig = scan_for_transits(time,Flux,period_guess-0.001,period_guess+0.00101,0.00001)
    time_filtered = t[np.isfinite(d)]
    depth_filtered = d[np.isfinite(d)]
    sigma_filtered = sig[np.isfinite(d)]
    
    return time_filtered[np.argmin(depth_filtered)]
    
    
