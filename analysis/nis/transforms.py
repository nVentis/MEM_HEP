from typing import Optional
from normflows.flows import Flow
import torch

class ScaleShift(Flow):
    def __init__(self, scale:torch.Tensor, shift:torch.Tensor,
                 limit_lower:Optional[torch.Tensor]=None, limit_upper:Optional[torch.Tensor]=None):
        super().__init__()
        
        self.scale = scale
        self.shift = shift
    
    def forward(self, z, context=None):
        z = z * self.scale + self.shift
        #log_det = torch.zeros(len(z), device=z.device)
        log_det = self.scale.log().sum()
        
        return z, log_det

    def inverse(self, z, context=None):
        z = z - self.shift
        z = z/self.scale
        #log_det = torch.zeros(len(z), device=z.device)
        log_det = -self.scale.log().sum()
        
        return z, log_det