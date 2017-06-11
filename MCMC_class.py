import numpy as np
import emcee
import kplr
import batman
import corner
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

def unwrap_self(arg, *emcee, **kwarg):
    return emcee[0].loglikelihood(arg, **kwarg)

class MCMC(object):

    def __init__(self, nwalkers, ndim,t, f, fe, bounds, numcores = 15):
        self.nwalkers = nwalkers
        self.ndim = ndim
        self.time = t
        self.flux = f
        self.error = fe
        self.bounds = bounds
        self.threads = numcores
    


    def loglikelihood(self, theta):
        """
        Compute the chi2 for a given set of parameters.
        """
        prior = self.lnPriorBounds(theta)
        chi2 = np.sum((self.flux-self.TransitModel(self.time, theta))**2/self.error**2)
        return prior - chi2/2
    def lnPriorBounds(self, theta):
        '''
        Set the priors for the MCMC for a given set of 
        paramters theta.
        theta = parameters
        bounds = list of tuples (lower bound, upper bound)

        returns 0 if theta in bounds and -inf if not in bound
        '''
        if np.any(theta<self.bounds[:,0]) | np.any(thet>self.bounds[:,0]):
            return -np.inf
        else:
            return 0
        

    def performMCMC(self, theta0,numIt, numPrune = 50):
        sampler = emcee.EnsembleSampler(self.nwalkers, self.ndim,
                                        unwrap_self, args = [self], threads = self.threads)

        p0 = [theta0 + 1e-4*np.random.randn(self.ndim) for i in range(self.nwalkers)]

        for i, result in enumerate(sampler.sample(p0, iterations=numIt)):
            if i == 0:
                print 'Starting MCMC'
            if (i+1) % 100 == 0:
                print("{0:5.1%}".format(float(i) / numIt))
        samples = sampler.chain[:, numPrune:, :].reshape((-1, self.ndim))
        print("Mean acceptance fraction: {0:.3f}"
                .format(np.mean(sampler.acceptance_fraction)))
        return samples

    '''Return results of the MCMC with given lower and upper bounds'''
    def results(self, samples, lower = 16, upper = 84):
        r = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
                             zip(*np.percentile(samples, [lower, 50, upper],
                                                axis=0)))
        params = []
        upper = []
        lower = []
        for result in r:
            params.append(result[0])
            upper.append(result[1])
            lower.append(result[2])
        return array(params), array(upper), array(lower)

    def cornerGraph(self, samples, label = ['Offset', 'Period', 'Radius', 'a', 'inc', 'e', 'peri', 'u1', 'u2']):
        fig = corner.corner(samples, bins=50,plot_datapoints=False,levels=[0.68,0.95],fill_contours=True,max_n_ticks=3,labels=label)
        plt.show()

    def plotTrans(self, params, width = 1):
        period = params[1]
        offset = params[0]

        timeFold = self.time%period
        fluxI = interp1d(self.time, self.flux, fill_value = 0)
        plt.plot(timeFold, fluxI(timeFold), 'b,')
        start = offset-width
        end = offset+width
        tTrans = np.linspace(start, end, 100)
        plt.plot(tTrans, self.TransitModel(tTrans, params), 'r')
        plt.xlim(start, end)
        plt.show()
        return
        
    def TracePlot(self):
        """
        Quick trace plot of the chain samples"""


    def TransitModel(self, time, parameters):
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