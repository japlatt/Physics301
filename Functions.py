import numpy as np
import emcee
import kplr
import transit

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