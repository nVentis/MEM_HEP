from analysis.mc.interfaces import ImportanceSamplingIntegrator
from analysis import tnp

class BinnedImportanceSamplingIntegrator(ImportanceSamplingIntegrator):
    def __init__(self, bins_per_dim:int=10, *args, **kwargs):
        super().__init__(*args, **kwargs)
        
        self.bins_per_dim = bins_per_dim
        bin_bounds = []
        
        for i in range(self.dims):
            boundary_arr = []
            c_cur = self.boundaries[i][0]
            c_step = (self.boundaries[i][1] - self.boundaries[i][0])/bins_per_dim
            
            for j in range(bins_per_dim):
                boundary_arr.append([c_cur, c_cur + c_step])
                c_cur += c_step
            
            bin_bounds.append(boundary_arr)
            
        #print("boundaries", self.boundaries.shape, self.boundaries)
    
        self.bin_indices = tnp.arange(bins_per_dim)
        self.bin_bounds = tnp.array(tnp.array(bin_bounds))
        
        
    # Generate
    def importance_per_dim(self, dim:int):
        bins = self.bin_bounds[dim][:].T
        diffs = bins[1] - bins[0]
        importance_per_dim = ((1/self.bins_per_dim))*1/diffs # (self.boundaries[dim][1] - self.boundaries[dim][0])/diffs
        
        return importance_per_dim
    
    def sample(self, n_samples:int=1, with_importance:bool=False):
        samples = []
        importance = []
        
        for i in range(self.dims):
            bin_idx = tnp.random.choice(self.bin_indices, size=n_samples)
            bins = self.bin_bounds[i][bin_idx].T
            
            samples_per_dim = tnp.random.uniform(low=bins[0], high=bins[1])
            samples.append(samples_per_dim)
            
            if with_importance:
                diffs = bins[1] - bins[0]
                importance_per_dim = ((1/self.bins_per_dim))*1/diffs # (self.boundaries[i][1] - self.boundaries[i][0])/diffs
                importance.append(importance_per_dim)
                
        if with_importance:
            return tnp.stack(samples), tnp.prod(tnp.stack(importance), axis=0)
        else:
            return tnp.stack(samples)
    
    # Adapt
    def adapt(self, n_samples:int=1000000):
        self.adapted = True
        samples = self.sample(n_samples)
        results = self.integrand(samples)
        
        for j in range(self.dims):
            # surface per dim and bin
            k = 0 #bin_lower = 0

            sort_mask = tnp.argsort(samples[j])
            samples_per_dim = samples[j][sort_mask]
            results_per_dim = results[sort_mask]
            bounds_per_dim = self.boundaries[j]
            
            surf = tnp.trapz(results_per_dim, x=samples_per_dim)
            surf_per_bin = surf/self.bins_per_dim
            
            bin_vals = []
            
            for i in range(self.bins_per_dim):
                running_sum = 0
                k0 = k
                while running_sum < surf_per_bin and k < n_samples - 1:
                    running_sum += results_per_dim[k]*(
                        (samples_per_dim[k+1] if k > 0 else bounds_per_dim[1]) - 
                        (samples_per_dim[k] if k <= n_samples-1 else bounds_per_dim[0]) )
                    k += 1
                
                bin_vals.append(tnp.array([
                    self.boundaries[j][0] if i == 0 else samples_per_dim[k0],
                    self.boundaries[j][1] if i == self.bins_per_dim-1 else samples_per_dim[k]]))
                
                #print(f"{k0}-{k} : {bin_vals[-1][0]}-{bin_vals[-1][1]}")
                
            self.bin_bounds[j] = tnp.stack(bin_vals)