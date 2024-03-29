import matplotlib.pyplot as plt
from analysis.fit_funcs import fit_funcs
from analysis.import_data import import_data, filter_data
from analysis.plot_matplotlib import plot_hist
from math import sqrt, pow, pi, exp
from typing import Optional, List
from os import path as osp
from os import makedirs, remove
import numpy as np

def plot_transfer(data, name, plot_save_dir:Optional[str] = None, fit:Optional[str]="gauss", true_label = "parton", reco_label = "jet",
                  quantity="E", xlabel=r"$ΔE$ [GeV]", xlim=(-100,100), ylim=(0, 0.18), n_bins=128, binrange=None, fit_init=None,
                  suptitle=None, fit_skip=False, yscale="linear", single_title=None, plot_args={}, titles:Optional[List[str]]=None):
    from scipy.optimize import curve_fit,minimize
    
    popts = []
    figures = []
    
    fig, axes = plt.subplots(1, len(data), figsize=(6*len(data),8))
    fig.suptitle((name + r": $" + quantity + r"_{" + reco_label + r"}-" + quantity +  r"_{" + true_label + r"}$") if suptitle is None else suptitle, fontsize=18)
    
    if len(data) == 1:
        axes = [axes]
    
    for i in range(1,1+len(data)):
        df = data[i-1]
        
        def dbgauss_construct(E_p, E_j):
            def dbgauss(a1, a2, a4, a5, b1, b2, b3, b4, b5):
                p1 = a1 + E_p*b1
                p2 = a2 + E_p*b2
                p3 =      E_p*b3
                p4 = a4 + E_p*b4
                p5 = a5 + E_p*b5
                
                return (1/(sqrt(2*pi)*(p2+p3*p5)))*(np.exp(-((E_j-E_p)-p1)**2/(2*p2**2)) + p3*np.exp(-((E_j-E_p) - p4)**2/(2*p5**2)))
            
            return dbgauss
        
        def dbgauss_likelihood(E_p, E_j):
            dbgauss = dbgauss_construct(E_p, E_j)
            def likelihood(a1, a2, a4, a5, b1, b2, b3, b4, b5):
                return np.sum(np.log(dbgauss(a1, a2, a4, a5, b1, b2, b3, b4, b5)))
        
        popt = None
        pcov = None
        df = df[0] - df[1]
        
        if fit == "dbgauss":
            raise Exception("Not implemented")
            # Unbinned likelihood fit
            #minimize(dbgauss_likelihood())
        elif isinstance(fit, str):            
            fit_inits = {
                "laplace": [np.average(df), 5],
                "lorentz": [np.average(df), 2],
                "uniform": [0.1],
                "gauss": [np.average(df), np.std(df)],
                "dbgauss_old": [0, np.std(df), 0.2, 0, 2*np.std(df)],
                "log_normal": [np.log(np.mean(df)) - sqrt(np.std(df))/2, np.std(df)]
            }
            
            fit_func = fit_funcs[fit]
            fit_init_c = fit_init if isinstance(fit_init, list) else (fit_inits[fit] if fit_init is None else fit_init[i-1])
            
            y, bins = np.histogram(df, bins=n_bins, density=True, range=binrange)
            x = (bins[:-1] + bins[1:]) / 2
            
            if fit_skip == True:
                popt = fit_init_c
            else:
                popt, pcov = curve_fit(fit_func, x, y, p0=fit_init_c, maxfev = 8000)
        
            print(popt)
        
        popts.append(popt)
        args = { 'x': (f"{reco_label.title()} {i}") if titles is None else titles[i-1], 'fit_func': (lambda x: fit_func(x, *popt)) if fit is not None else None,
                'fit_opts': popt, 'bins': n_bins, 'xlim': xlim, 'ylim': ylim, 'ax': axes[i-1],
                 'xlabel': xlabel, 'title':(None if single_title is None else single_title), 'normalize':True, 'yscale': yscale, 'text_spacing_y': 0.15 }
        fig = plot_hist(df, **{ **args, **plot_args})
        fig.tight_layout()
        figures.append(fig)
        
    if len(popts) == 4 and (not 0 in [(1 if a is not None else 0) for a in popts]):
        print("SEPTF: A(1+2) : B(3+4)")
        A = np.array([
            popts[0],
            popts[1]
        ])

        B = np.array([
            popts[2],
            popts[3]
        ])

        print(A.mean(axis=0))
        print(B.mean(axis=0))
    
    return figures, popts

def plot_transfer_from_df(data, name = "", plot_save_dir:Optional[str] = None, fit = "gauss", yscale="linear"):
    df = []
    for i in range(1, 5):
        df.append((data[f"jet{i}_e"], data[f"parton{i}_e"]))
    
    plot_transfer(df, name, plot_save_dir=plot_save_dir, fit=fit, yscale=yscale)