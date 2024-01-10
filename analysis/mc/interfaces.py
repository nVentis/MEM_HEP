from typing import Callable, Iterable, List

class IntegratorInterface:
    def __init__(self, integrand:Callable, dims:int, boundaries:List[List[float]], mode:int=0, **kwargs) -> None:
        self.integrand = integrand
        self.dims = dims
        self.boundaries = boundaries
        self.mode = mode
    
    def integrate(self, **kwargs):
        raise NotImplementedError()
    
class ImportanceSamplingIntegrator(IntegratorInterface):
    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self.adapted = False
        
    def adapt(self, **kwargs):
        raise NotImplementedError()
    
    def sample(self, n_samples:int=1, with_importance:bool=False, **kwargs):
        raise NotImplementedError()