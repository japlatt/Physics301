import numpy as np
import emcee
import kplr
import batman

def Fw(t,t0,F0,sigma):
    """
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


def correlate(flux, n_coors = 80, std_filter = 5):
    coor = np.zeros(n_coors)
    mean = np.nanmean(flux)
    standard_dev = np.nanstd(flux)
    for delta in xrange(n_coors):
        I1 = np.copy(flux)
        I2 = array(list(flux[delta:]) + list(flux[0:delta]))

        inds1 = (I1 >= mean - std_filter*standard_dev) & \
                (I1 <=  mean + std_filter*standard_dev) & \
                 (I1 >= 0.5)
        inds2 = (I2 >= mean - std_filter*standard_dev) & \
                (I2 <=  mean + std_filter*standard_dev) & \
                (I2 >= 0.5)
        inds = inds1 & inds2
        N = sum(inds)
        coor[delta] = sum((I1[inds]*I2[inds]))/N
    return coor
    
def Cn(time,flux,deltat):
    """
    Compute the correlation function for a given delta T scanned over all times
    """
    C = 0.0
    N = 0
    for i in range(len(time)):
        Fi = flux[i]
        Fipdt = Fw(time[i]+deltat,time,flux,1./24.)
        if not np.isclose(Fipdt,0):
            C += Fi*Fipdt
            N += 1
    
    C /= float(N)
    return C

