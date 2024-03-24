import torch
import numpy as np
from typing import Optional,Union
from normflows.distributions import BaseDistribution

class Uniform(BaseDistribution):
    """
    Multivariate uniform distribution
    """

    def __init__(self, shape, low=-1.0, high=1.0, device:Union[torch.device,str]='cpu', deterministic:bool=False):
        """Uniform distributions with the same bounds along all dimensions given by shape

        Args:
            shape (_type_): _description_
            low (float, optional): _description_. Defaults to -1.0.
            high (float, optional): _description_. Defaults to 1.0.
            device (Union[torch.device,str], optional): _description_. Defaults to 'cpu'.
            deterministic (bool, optional): if True, will return deterministic samples along uniform grid. Defaults to False.
        """
        super().__init__()
        if isinstance(shape, int):
            shape = (shape,)
        if isinstance(shape, list):
            shape = tuple(shape)
        self.shape = shape
        self.d = int(np.prod(shape))
        self.register_buffer('log_prob_val', tensor=torch.tensor(-self.d * np.log(high - low)))
        self.register_buffer('low', tensor=torch.tensor(low))
        self.register_buffer('high', tensor=torch.tensor(high))
        self.device = device
        self.deterministic = deterministic

    def forward(self, num_samples=1, context=None):
        if not self.deterministic:
            eps = torch.rand(
                (num_samples,) + self.shape, dtype=self.low.dtype, device=self.device
            )
            z = self.low + (self.high - self.low) * eps
        else:           
            z = (torch.linspace(self.low, self.high, steps=(num_samples+1), dtype=self.low.dtype, device=self.low.device)[:-1])
            z = z[:, None].expand(num_samples, self.d)
            
        log_p = self.log_prob_val * torch.ones(num_samples, device=self.device)
        return z, log_p

    def log_prob(self, z, context=None):
        log_p = self.log_prob_val * torch.ones(z.shape[0], device=z.device)
        out_range = torch.logical_or(z < self.low, z > self.high)
        ind_inf = torch.any(torch.reshape(out_range, (z.shape[0], -1)), dim=-1)
        if len(ind_inf):
            print('ind_inf', ind_inf)
        log_p[ind_inf] = -np.inf
        return log_p

class HyperUniform(BaseDistribution):
    """
    Uniform distribution along customizable bounds per dimensions
    """

    def __init__(self, low, high, dtype=torch.float, device:Optional[Union[str,torch.device]]=None):
        """Constructor

        Args:
          low: Lower bound of uniform distribution
          high: Upper bound of uniform distribution
        """
        
        low = torch.tensor(low).to(dtype)
        high = torch.tensor(high).to(dtype)
        
        if device is not None:
            low = low.to(device)
            high = high.to(device)
        
        assert(low.shape == high.shape)
        
        shape = low.shape[0]
        
        super().__init__()
        if isinstance(shape, int):
            shape = (shape,)
        if isinstance(shape, list):
            shape = tuple(shape)
            
        self.shape = shape
        self.d = np.prod(shape)
        self.low = low
        self.high = high
        self.log_prob_val = -(self.high - self.low).log().sum()

    def forward(self, num_samples=1, context=None):
        eps = torch.rand(
            (num_samples,) + self.shape, dtype=self.low.dtype, device=self.low.device
        )
        z = self.low + (self.high - self.low) * eps
        log_p = self.log_prob_val * torch.ones(num_samples, device=self.low.device)
        return z, log_p

    def log_prob(self, z, context=None):
        return torch.ones(z.shape[0], device=z.device)
        log_p = self.log_prob_val * torch.ones(z.shape[0], device=z.device)
        out_range = torch.logical_or(z < self.low, z > self.high)
        ind_inf = torch.any(torch.reshape(out_range, (z.shape[0], -1)), dim=-1)
        log_p[ind_inf] = -np.inf
        return log_p