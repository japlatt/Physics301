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
        I2 = array(list(flux[delta:]) + list(flux[0:delta]))

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

