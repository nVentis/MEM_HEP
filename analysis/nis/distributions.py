import torch
import numpy as np
from normflows.distributions import BaseDistribution

class HyperUniform(BaseDistribution):
    """
    Multivariate uniform distribution
    """

    def __init__(self, low, high, dtype=torch.float):
        """Constructor

        Args:
          low: Lower bound of uniform distribution
          high: Upper bound of uniform distribution
        """
        
        low = torch.tensor(low).to(dtype)
        high = torch.tensor(high).to(dtype)
        
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
        self.log_prob_val = -torch.sum(torch.log(self.high - self.low))

    def forward(self, num_samples=1, context=None):
        eps = torch.rand(
            (num_samples,) + self.shape, dtype=self.low.dtype, device=self.low.device
        )
        z = self.low + (self.high - self.low) * eps
        log_p = self.log_prob_val * torch.ones(num_samples, device=self.low.device)
        return z, log_p

    def log_prob(self, z, context=None):
        log_p = self.log_prob_val * torch.ones(z.shape[0], device=z.device)
        out_range = torch.logical_or(z < self.low, z > self.high)
        ind_inf = torch.any(torch.reshape(out_range, (z.shape[0], -1)), dim=-1)
        log_p[ind_inf] = -np.inf
        return log_p