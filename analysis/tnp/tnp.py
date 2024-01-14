from types import SimpleNamespace
from typing import Optional
import numpy as np
import torch

conv_dict = {
    # numpy mode
    0: {
        'arange': np.arange,
        'array': np.array,
        'prod': np.prod,
        'trapz': np.trapz,
        'argsort': np.argsort,
        'sum': np.sum,
        'exp': np.exp,
        'stack': np.stack,
        'ones': np.ones,
        'zeros': np.zeros,
        'random': {
            'choice': np.random.choice,
            'uniform': np.random.uniform,
        }
        
    },
    # torch mode
    1: {
        'arange': torch.arange,
        'array': torch.tensor,
        'prod': lambda tensor, axis=0: torch.prod(tensor, dim=axis),
        'trapz': torch.trapezoid,
        'argsort': torch.argsort,
        'sum': torch.sum,
        'exp': torch.exp,
        'stack': torch.stack,
        'ones': torch.ones,
        'zeros': torch.zeros,
        'random': {
            'choice': lambda vals, size: torch.multinomial(torch.ones(len(vals))*1/len(vals), num_samples=size, replacement=True),
            'uniform': lambda low, high: (high-low)*torch.rand(len(low)) + low,
        }
    }
}

__cur_mode__ = 0
tnp = SimpleNamespace(**conv_dict[__cur_mode__])

def tnp_mode(mode:Optional[int]=None):
    global __cur_mode__
    
    if mode is None:
        return __cur_mode__
    else:
        __cur_mode__ = mode
        
        if not (mode in conv_dict):
            raise NotImplemented(f'Mode {mode} not supported')
        
        def merge(base, items):
            for k, v in items:
                if isinstance(v, dict):
                    setattr(base, k, SimpleNamespace(**v))
                    #merge(getattr(base, k), v)
                else:
                    setattr(base, k, v)
        
        merge(tnp, conv_dict[mode].items())
        
        return tnp

tnp_mode(__cur_mode__)