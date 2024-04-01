import torch
import numpy as np
import normflows as nf
import pyro.distributions as dist
import pandas as pd
import matplotlib.pyplot as plt
from analysis.plot_matplotlib import plot_hist
from typing import Optional, Callable

def plot_proposal(nfm, n_samples:int=1024):
    with torch.no_grad():
        if isinstance(nfm, nf.NormalizingFlow):
            samples = nfm.sample(n_samples)[0]
        else:
            samples = nfm.sample(n_samples)
            
        samples = samples.cpu().numpy()

        all_max = np.max(samples).max()
        all_min = np.min(samples).min()
        print(f'Using inferred bounds {all_min} to {all_max}')
        
        #print((all_min, all_max))
        dfbase = {}
        for i in range(len(samples.T)):
            dfbase[f'x{i}'] = samples.T[i]

        plot_hist(pd.DataFrame(dfbase), normalize=True, xlim=(all_min, all_max), xlim_binning=(all_min, all_max))

def plot_integrand(nfm, func, n_samples:int=2048, vwindow=None, y_log:bool=False):
    with torch.no_grad():    
        if isinstance(nfm, nf.NormalizingFlow):
            samples = nfm.sample(n_samples)[0]
        else:
            samples = nfm.sample(n_samples)
        
        samples = samples.cpu()
        results = func(samples).cpu().numpy()
        
        for i in range(len(samples.T)):
            fig, ax = plt.subplots()
            #ax.hist(samples.T[i], alpha=0.4, label=f"Dim {i}", bins=64)
            ax.scatter(x=samples.T[i].numpy(),
                    y=results)
            
            if vwindow is not None:
                if len(vwindow) > 2 and len(vwindow) == len(samples.T):
                    ax.set_ylim(vwindow[i])
                else:
                    ax.set_ylim(vwindow)
                
            ax.set_title('Integrand projection')
            ax.set_xlabel(f"Dim {i}")
            ax.set_ylabel('Integrand value')
            plt.show()
        
        fig, ax = plt.subplots()
        ax.hist(results, bins=32)
        ax.set_xlabel('Integrand value')
        if y_log:
            ax.set_yscale('log')
        ax.set_title('All results')
    
def get_result(nfm, func:Callable, n_samples:int=1024):
    with torch.no_grad():
        if isinstance(nfm, nf.NormalizingFlow):
            samples, log_p = nfm.sample(n_samples)
            test = log_p.exp()#.cpu().numpy()   
        else:
            samples = nfm.sample(n_samples)
            test = nfm.log_prob_val.exp()
        
        true = func(samples)
        mean, var = torch.mean(true/test).item(), torch.var(true/test).item()
        
    return mean, var