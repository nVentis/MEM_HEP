import torch
import numpy as np

def gaussian(x, alpha:float=.2, center=.5):
    prefac = 1./(alpha*np.sqrt(np.pi))**x.shape[1]
    exp1 = -1.*torch.sum(((x-center)**2)/alpha**2, dim=-1)
    
    return prefac*exp1.exp()

def camel(x:torch.Tensor, alpha:float=.2, centers=[1./3., 2./3.]):
    """Two Gaussians centered at (1/3, 2/3) in each dimension.

    Integral equals
        (0.5*(erf(1/(3*alpha)) + erf(2/(3*alpha)) ))** ndims
        
    Args:
        x (torch.Tensor): _description_
        alpha (float, optional): _description_. Defaults to 0.2.

    Returns:
        _type_: _description_
    """
    
    prefac = .5/(alpha*np.sqrt(np.pi))**x.shape[1]
    exp1 = -1.*(((x-centers[0])/alpha)**2).sum(axis=1)
    exp2 = -1.*(((x-centers[1])/alpha)**2).sum(axis=1)
    
    if isinstance(x, torch.Tensor):
        return prefac*(exp1.exp()+exp2.exp())
    else:
        return prefac*(np.exp(exp1)+np.exp(exp2))