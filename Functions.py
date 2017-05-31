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
    lcs = planet.get_lightcurves(short_cadence=False)
    time = np.zeros(0)
    flux = np.zeros(0)
    ferr = np.zeros(0)
    
    for lc in lcs:
        with lc.open() as f:
            hdu_data = f[1].data()
            time = np.append(time,hdu_data["time"])
            flux = np.append(flux,hdu_data["PDCSAP_FLUX"])
            ferr = np.append(ferr,hud_data["PDCSAP_FLUX_ERR"])
            
    Time_cleaned = Time[np.isfinite(flux)]
    Ferr_cleaned = ferr[np.isfinite(flux)]
    Flux_cleaned = flux[np.isfinite(flux)]
    
    Ferr_norm = Ferr_cleaned / np.median(Flux_cleaned)
    Flux_norm = Flux_cleaned / np.median(FLux_cleaned)
    
    return Time_cleaned , Flux_norm , Ferr_norm


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
    for j,t in enumerate(deltaT):
        if np.isclose(t%0.25,0): print t
        for i in range(len(bins)-1):
            inds = np.logical_and(((time%t)/t>bins[i]),(time%t)/t<=bins[i+1])

            Flux[i] = np.median(flux[inds])
        depth[j] = np.min(Flux)
        
    return depth ,deltaT
    
def Identify_transits(deltaT,depth):
    """
    parse the full array of deltaT and depth, and identify the individual transit candidates, including
    their time of peak signal and the depth of that signal
    
    return cleaned list of transit candidates (any depth greater than 1 standard deviation)
    """
    
    i = 10
    
    depth_processed = (depth-np.median(depth))/np.std(depth)
    transit_list = []
    print np.sum(depth_processed <= -1000)
    
    while i<len(deltaT):
        if i%1000 ==0: print depth_processed[i]
        if depth_processed[i] <= -1:
            # transit candidate, search for brightest nearby signal
            Transit_time = 0.0
            Transit_depth = np.inf
            for j in range(i-10,i):
                if depth[j] < Transit_depth:
                    Transit_time = deltaT[j]
                    Transit_depth = depth[j]
            # have scanned from behind the first detected location, now want to scan ahead
            # do so until the processed signal is back above 1 sigma
            j = 0
            while depth_processed[j+i] <= -1:
                if depth[j+i] < Transit_depth:
                    print Transit_depth,depth[j+i]
                    Transit_time = deltaT[j+i]
                    Transit_depth = depth[j+i]
                j += 1
            transit_list.append([Transit_time,Transit_depth])
            i += j
            
        else:
            i +=1
    
    Transit_array = np.zeros([len(transit_list),2])
    for i in range(len(transit_list)):
        Transit_array[i,0] = transit_list[i][0]
        Transit_array[i,1] = transit_list[i][1]
    return Transit_array
        
    
def Get_scan_info(deltaT,depth,firstN,Flux,Time):
    """
    For the firstN strongest transits (as measured by depth), calculate the following 
    two quantities:
    
    Ntransits:          The number of transits that occur when the lightcurve is folded on this timescale
                        This is useful for discriminating real transits from harmonics.
    
    Transitdepthstd:    The standard deviation of the measured fluxes at the peak of the transit signal.  
                        We expect this should be larger if the timescale deltaT is a fraction of the true
                        transit signal, because it will place transits in resonance, but also non-transits
                        at the same location.
    """
    depthsorted = [y for (y,x) in sorted(zip(depth,deltaT))]
    Timesorted  = [x for (y,x) in sorted(zip(depth,deltaT))]
    
    data = np.zeros([firstN,3])
    
    for j in range(firstN):
        t = Fsorted[i]
        Ncycles = 0
        for i in range(len(bins)-1):
            inds = np.logical_and(((Time[:50000]%t)/t>bins[i]),(Time[:50000]%t)/t<=bins[i+1])

            Flux[i] = np.mean(Flux[:50000][inds])
            Fluxstd[i] = np.std(Flux[:50000][inds])

        Df = Flux[:-1] - Flux[1:]
        for i in range(1,len(Df)):
            if (Df[i]>3*np.std(Df)) & (Df[i-1]<3*np.std(Df)):
                Ncycles += 1
            
        #data[j,0] = 
    


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

