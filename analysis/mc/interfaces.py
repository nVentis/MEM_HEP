from typing import Callable, Iterable

class IntegratorInterface:
    def __init__(self, integrand:Callable, boundaries:list[list[float]], **kwargs) -> None:
        self.integrand = integrand
        self.dims = len(boundaries)
        self.boundaries = boundaries
    
    def integrate(self, **kwargs):
        raise NotImplementedError()
    
class ImportanceSamplingIntegrator(IntegratorInterface):
    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self.adapted = False
        
    def adapt(self, *args, **kwargs):
        raise NotImplementedError()
    
    def sample(self, n_samples:int=1, with_importance:bool=False, *args, **kwargs):
        raise NotImplementedError()
    
    # Integrate
    def integrate(self, n_samples:int=10000):
        samples, importance = self.sample(n_samples, with_importance=True)
        results = self.integrand(samples)
        
        V = 1
        for j in range(self.dims):
            V = V*(self.boundaries[j][1] - self.boundaries[j][0])
                
        if False:#not self.adapted:
            res = results.sum()/n_samples
                
            return res*V
        else:
            f = results
            #p = 1/V
            #q = importance/V
            p = importance
            
            #print("p", p)
            #print("p/q", p/q)
            #print("f", f)
            #print("q", q)
            
            #return (f*p/q).sum().sum()/n_samples
            return (f/p).sum().sum()/n_samples