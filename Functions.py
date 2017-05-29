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

def correlate(flux, n_coors = 80, std_filter = 5):
    coor = np.zeros(n_coors)
    mean = np.nanmean(flux)
    standard_dev = np.nanstd(flux)
    for delta in xrange(n_coors):
        I1 = copy(flux)
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





