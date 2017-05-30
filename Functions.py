import numpy as np
import emcee
import kplr
from scipy.optimize import interp1d
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
        self.VECTORS = np.random.random([len(priorup),N])
        for i in range(len(ub)):
                self.VECTORS[i,:] *= (ub[i]-lb[i])
                self.VECTORS[i,:] += lb[i]
                
        self.DONORS = np.zeros(self.VECTORS.shape)
        self.TRIALS = np.zeros(self.VECTORS.shape)
        self.CURRENT_FX = np.zeros(self.VECTORS.shape[1])
        self.StrategyProbs = np.array([0.25,0.25,0.25,0.25])
        return
        
    def Mutation(self):
        """
        Mutate the vectors.
        """
        for i in range(self.VECTORS.shape[1]):
            options = range(self.VECTORS.shape[1])
            options.remove[i]
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
        for i in range(vectors.shape[1]):
            f1 = self.CURRENT_FX[i]
            f2 = ftominimize(x,y,yerr,trials[:,i])
            if f2 <= f1:
                self.VECTORS[:,i] = self.TRIALS[:,i]
            else:
                pass
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
        for i in range(Niter):
            self.Mutation()
            self.Recombination()
            self.Selection()
            
            if np.mean(np.std(self.VECTORS,axis=1)) < convThresh:
                break
                
        print "Optimization has Converged"
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

