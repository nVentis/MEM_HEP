import numpy as np
from typing import List, Union, Tuple

# From iflow
def variance_weighted_result(means:Union[List,np.ndarray], stddevs:Union[List,np.ndarray])->Tuple[float, float]:
    """ Computes weighted mean and stddev of given means and
        stddevs arrays, using Inverse-variance weighting
    """
    
    if isinstance(means, list): means = np.array(means)
    if isinstance(stddevs, list): stddevs = np.array(stddevs)
    
    assert np.size(means) == np.size(stddevs)
    assert means.shape == stddevs.shape
    variance = 1./np.sum(1./stddevs**2, axis=-1)
    mean = np.sum(means/(stddevs**2), axis=-1)
    mean *= variance
    return mean, np.sqrt(variance)
