import numpy as np
import emcee
import kplr
#import batman

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

<<<<<<< HEAD

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
        Set the initalization region
        """
        return
    def Mutation(self):
        """
        Mutate the vectors, using the 4 rules.
        """
        return
    def Recombination(self):
        """
        Recombine the donor vectors with the current vectors
        to create the trial vectors
        """
        return
    def Selection(self,likelihoodfunction):
        """
        call likelihoodfunction with the trial vectors, accept
        if it improves the fit, and replace existing vectors
        with the trial vectors.
        """
        return
    def UpdateStrategyProbs(self):
        """
        Update the probabilities that a given strategy will be attempted, as well
        as CR and F for recombination.
        """
        return
    def Optimize(Niter,convThresh):
        """
        Run the optimization for either Niter iterations, or until the rms scatter 
        over all the dimensions is lower than convThresh
        """
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
        I1 = copy(flux)
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

