import numpy as np
import emcee
import kplr
import batman
import corner
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

def unwrap_self(arg, *emcee, **kwarg):
    return emcee[0].lnPost(arg, **kwarg)

class MCMC(object):

    def __init__(self, nwalkers, ndim,t, f, fe, bounds, numcores = 20):
        self.nwalkers = nwalkers
        self.ndim = ndim
        self.time = t - min(t)
        self.flux = f
        self.error = fe
        self.bounds = bounds
        self.threads = numcores
    

    def lnPost(self, theta):
        prior = self.lnPriorBounds(theta)
        if ~np.isfinite(prior):
            return prior
        
        return prior + self.loglikelihood(theta)

    def loglikelihood(self, theta):
        """
        Compute the chi2 for a given set of parameters.
        """
        chi2 = np.sum(((self.flux - self.TransitModel68C(self.time, theta))**2)/self.error**2)
        return -chi2/2
            
    def lnPriorBounds(self, theta):
        '''
        Set the priors for the MCMC for a given set of 
        paramters theta.
        theta = parameters
        bounds = list of tuples (lower bound, upper bound)

        returns 0 if theta in bounds and -inf if not in bound
        '''
        if np.any(theta < self.bounds[:,0]) or np.any(theta > self.bounds[:,1]):
            return -np.inf
        else:
            return 0

    #def addSamples(chain):





    def performMCMC(self, theta0, numIt, name = None):
        sampler = emcee.EnsembleSampler(self.nwalkers, self.ndim,
                                        unwrap_self, args = [self], threads = self.threads)
        #self.lc = self.initTransitModel(self.time, theta0)

        p0 = [theta0 + 1e-4*np.random.randn(self.ndim) for i in range(self.nwalkers)]

        for i, result in enumerate(sampler.sample(p0, iterations=numIt)):
            if i == 0:
                print 'Starting MCMC'
            if (i+1) % 10 == 0:
                print("{0:5.1%}".format(float(i) / numIt))
                print("Mean acceptance fraction: {0:.3f}".format(np.mean(sampler.acceptance_fraction)))
            if (i+1) % 100 == 0:
                if name != None:
                    np.save(name, sampler.chain)
                


        print("Mean acceptance fraction: {0:.3f}"
                .format(np.mean(sampler.acceptance_fraction)))
        self.samples = sampler.chain.reshape((-1, self.ndim))
        self.origSample = sampler.chain
        self.sampler = sampler

    def pruneSamples(self, numPrune):
        #self.samples =(self.samples.reshape(self.nwalkers, numPrune, self.ndim)[:,numPrune:, :]).reshape(-1, self.ndim)
        self.samples = self.origSample[:, numPrune:, :].reshape((-1, self.ndim))

    def chain(self):
        return self.samples

    def origChain(self):
        return self.origSample

    def saveSamples(self, name):
        np.save(name, self.samples)

    def TransitModel68(self, time, parameters):
        params     = batman.TransitParams()
        params.t0  = parameters[0]                      # time of inferior conjunction
        params.per = parameters[1]                      # orbital period
        params.rp  = parameters[2]                      # planet radius (in units of stellar radii)
        params.a   = parameters[3]                      # semi-major axis (in units of stellar radii)
        params.inc = 87.68                              # orbital inclination (in degrees)
        params.ecc = 0.02                               # eccentricity
        params.w   = 90                                 # longitude of periastron (in degrees)
        params.u   = [0.3908288, 0.263158]              # limb darkening coefficients [u1, u2]
        params.limb_dark = "quadratic"                  # limb darkening model
        
        model = batman.TransitModel(params, time)
        return model.light_curve(params)
    def TransitModel68C(self, time, parameters):
        params     = batman.TransitParams()
        params.t0  = parameters[0]                      # time of inferior conjunction
        params.per = parameters[1]                      # orbital period
        params.rp  = parameters[2]                      # planet radius (in units of stellar radii)
        params.a   = parameters[3]                      # semi-major axis (in units of stellar radii)
        params.inc = 86.93                              # orbital inclination (in degrees)
        params.ecc = 0.0                                # eccentricity
        params.w   = 90                                 # longitude of periastron (in degrees)
        params.u   = [0.3908288, 0.263158]              # limb darkening coefficients [u1, u2]
        params.limb_dark = "quadratic"                  # limb darkening model

        model = batman.TransitModel(params, time)
        return model.light_curve(params)

    '''Return results of the MCMC with given lower and upper bounds'''
    def results(self, lower = 16, upper = 84):
        samples = self.samples
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
        return np.array(params), np.arrIay(upper), np.array(lower)

    def cornerGraph(self, label = ['Offset', 'Period', 'Radius', 'a', 'inc', 'e', 'peri', 'u1', 'u2']):
        fig = corner.corner(self.samples, bins=50,
            plot_datapoints=False,levels=[0.68,0.95],
            fill_contours=True,max_n_ticks=3,labels=label,
            quantiles=[0.16, 0.5, 0.84], show_titles=True, title_kwargs={"fontsize": 15})
        return fig

    def plotTrans(self, params, width = 0.5):
        period = params[1]
        offset = params[0]

        fig = plt.figure(figsize = (10,7))
        timeFold = self.time%period
        fluxFold = self.flux%period
        start = offset-width
        end = offset+width
        tTrans = np.linspace(start, end, 100)
        #fluxI = interp1d(self.time, self.flux, fill_value = 0)
        plt.plot(timeFold, fluxFold, 'k,')

        plt.plot(tTrans, self.TransitModel(tTrans, params), 'r')
        plt.xlim(start, end)
        return fig


    def tracePlot(self):
        def plotLine(ax, samples, k, label = ['Offset', 'Period', 'Radius', 'a', 'inc', 'e', 'peri', 'u1', 'u2']):
            length = samples.shape[0]
            for i in xrange(length):
                ax.plot(samples[i, :, k], linewidth = 0.5, alpha = 0.3)

        ax.set_title(label[k])
        chain = self.origChain
        fig, ax = plt.subplots(3,3, figsize = (11,8))
        plotLine(ax[0,0], chain, 0)
        plotLine(ax[0,1], chain, 1)
        ax[0,1].set_ylim(5.39875, 5.39876)
        plotLine(ax[0,2], chain, 2)
        plotLine(ax[1,0], chain, 3)
        plotLine(ax[1,1], chain, 4)
        plotLine(ax[1,2], chain, 5)
        plotLine(ax[2,0], chain, 6)
        plotLine(ax[2,1], chain, 7)
        plotLine(ax[2,2], chain, 8)
        return fig

    def setParams(self, parameters):
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
        return params

    def TransitModel(self, time, parameters):
        """
        Use the Batman transit modeling code to produce
        a transit model brightness (normalized to 1), at
        each time point.
        """
        params = self.setParams(parameters)
        model = batman.TransitModel(params, time)
        return model.light_curve(params)

    def initTransitModel(self, time, parameters):
        """
        Use the Batman transit modeling code to produce
        a transit model brightness (normalized to 1), at
        each time point.
        """
        
        params = self.setParams(parameters)
        model = batman.TransitModel(params, time)
        return model
    def lightCurve(self, model, parameters):
        params = self.setParams(parameters)
        return model.light_curve(params)



