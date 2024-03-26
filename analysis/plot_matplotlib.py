import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
import matplotlib.pylab as pylab
import pandas as pd
from matplotlib import rcParams as rcp
from matplotlib.backends.backend_pdf import PdfPages
from typing import Optional, Union, Callable, Dict, List
try:
    from typing import Literal
except ImportError:
    from typing_extensions import Literal

from math import sqrt
#from scipy.stats import chisquare

rcp['hatch.linewidth'] = 0.5  # previous pdf hatch linewidth
rcp['font.family'] = 'monospace' # per default use monospace fonts

def format_st(f, t = 0.01):
    return f"{f:.2f}" if f >= t else f"<{t}"

settings = {
    'colorpalette': None
}

def get_colorpalette():
    if settings['colorpalette'] is None:
        #colorpalette = ["tab:blue", "tab:red", "y", "tab:pink", "tab:cyan", "tab:olive"]
        prop_cycle = plt.rcParams['axes.prop_cycle']
        colorpalette = prop_cycle.by_key()['color']
        settings['colorpalette'] = colorpalette

    return settings['colorpalette']
    
def set_colorpalette(colorpalette):
    settings['colorpalette'] = colorpalette


def fontsize(fs):
    pylab.rcParams.update({
        'legend.fontsize': fs,
        'axes.labelsize': fs,
        'axes.titlesize': fs,
        'xtick.labelsize': fs,
        'ytick.labelsize': fs})

def export_figures(filename, figs=None, dpi=200):
    pp = PdfPages(filename)
    if figs is None:
        figs = [plt.figure(n) for n in plt.get_fignums()]
        #print(len(figs))
        
    for fig in figs:
        fig.savefig(pp, format='pdf')
    pp.close()

def plot_hist(data:Union[dict,pd.DataFrame], x:Optional[Union[str,list]]=None,
              fit_func:Optional[Callable]=None, fit_opts:Optional[Dict] = None,
              labels:Optional[List[str]]=None, colorpalette=None, bins=128,
              same_bins:bool=True, xlim_binning:Optional[Union[list,tuple]]=None, xlim:Optional[Union[list,tuple]]=None, ylim=None,
              xlabel:Optional[str] = None, ylabel:Optional[str]=None,
              normalize=False, filter_nan:bool=False,
              title:Optional[str]=None, ax=None,
              text_start_x:float=0.965, text_start_y:float=0.97, text_spacing_y:float=0.22,
              xscale:Literal['linear', 'log']='linear', yscale:Literal['linear', 'log']="linear",
              fontsize:Optional[Union[str, int]]=14, legendsize = None, titlesize:Union[int, str]=15,
              ticksize_minor:int=10, ticksize_major=None,
              figsize:tuple=(8,6), figdpi:int=100, scientific_stats:bool=False):
    """_summary_
    
    text_spacing_y: 0.11 for high-res
    

    Args:
        data (Union[dict,pd.DataFrame]): Can be a pd.DataFrame or simple dict (preferrable for columns of unequal size)
        x (Optional[Union[str,list]]): one column or a list of columns in data to be histogrammed.
        fit_func (Optional[function], optional): only supported if one column is to be plotted, i.e. x is a string
        fit_opts (Optional[dict], opional): only used for printing
        labels (_type_, optional): _description_. Defaults to None.
        colorpalette (_type_, optional): _description_. Defaults to None.
        bins (int, optional): _description_. Defaults to 128.
        same_bins (bool): if True, the same bin ranges are forced for all plotted histograms. Defaults to True.
        xlim_binning (list, optional): _description_. Defaults to False.
        xlim (Optional[list], optional): _description_. Defaults to None.
        ylim (_type_, optional): _description_. Defaults to None.
        xlabel (Optional[str], optional): _description_. Defaults to None.
        ylabel (Optional[str], optional): _description_. Defaults to None.
        normalize (bool, optional): _description_. Defaults to False.
        title (Optional[str], optional): _description_. Defaults to None.
        ax (_type_, optional): _description_. Defaults to None.
        filter_nan (): Whether or not to filter out nan data; useful for DataFrames of unequal sizes. Defaults to False.
        text_start_x (float, optional): _description_. Defaults to 1.02.
        text_start_y (float, optional): _description_. Defaults to 1.05.
        text_spacing_y (float, optional): Height of each textbox . Defaults to 0.11.
        xscale (literal, optional): _description_. Defaults to "linear".
        yscale (literal, optional): _description_. Defaults to "linear".
        fontsize (Optional[str], optional): _description_. Defaults to None.
        ticksize_major (int, optional): ticksize_minor+2 if None. Defaults to None
    """
    
    if ax == None:
        fig, ax = plt.subplots()
        fig.set_dpi(figdpi)
        fig.set_figwidth(figsize[0])
        fig.set_figheight(figsize[1])
    else:
        fig = plt.gcf()
    
    if colorpalette is None:
        colorpalette = get_colorpalette()
    
    # Force conversion of dict to DataFrame
    if isinstance(data, dict):
        data = pd.DataFrame(data)
    elif isinstance(data, np.ndarray):
        data = pd.DataFrame({ (x if isinstance(x, str) else "Data"): data })
    
    x_in = x
    
    if x is None:
        x = list(data.keys())
    
    """
    if isinstance(data, dict):
        if len(data.keys()) == 1:
            columns = [None]
            data = data[data.keys()[0]]
            xlim_view = [0.98*data.min(), 1.02*data.max()] if xlim is None else xlim
        else:
            columns = x
            if xlim_view is None:
                xlim_view = [float('inf'), float('-inf')]
                for x in data:
                    xlim_view[0] = min(xlim_view[0], 0.98*data[x].min())
                    xlim_view[1] = max(xlim_view[1], 1.02*data[x].max())
                
    else: # DataFrame
    """
    
    xlim_view = xlim if xlim is not None else None
    
    
    if isinstance(x, str):
        x = [x]

    if len(list(data.shape)) == 1:
        columns = [None] # In this case, data is assumed to contain just one column of data, which is to be histogrammed
        if xlim_view is None:
            xlim_view = [0.98*data.min(), 1.02*data.max()]
    else:
        columns = [x] if isinstance(x, str) else x
        if xlim_view is None:
            xlim_view = [0.98*data[x].min().min(), 1.02*data[x].max().max()]
                
    # If same_bins=True, infer limits and impose xlim_binning
    if same_bins and xlim_binning is None:
        if isinstance(data, pd.DataFrame):
            xlim_binning = xlim_view
        else:
            xlim_binning = [np.min(np.min(data)), np.max(np.max(data))]

    for i in range(len(columns)):
        print("asd")
        column = columns[i]
        values = data if column is None else data[column]
        
        if filter_nan == True:
            values = values[~np.isnan(values)]
            
        # Limits
        min_val = xlim_binning[0] if (xlim_binning is not None) else np.min(values)
        max_val = xlim_binning[1] if (xlim_binning is not None) else np.max(values)
        
        bin_edges = np.linspace(min_val, max_val, bins + 1)
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
        
        # Filter for x_lim_binning
        stat_values = values
        if xlim_binning is not None:
            stat_values = stat_values[(stat_values >= xlim_binning[0]) & (stat_values <= xlim_binning[1])]
        
        # Additional behavior if only one column is to be plotted
        h_name  = (x if isinstance(x, str) else "Data") if column is None else (column if labels is None else labels[i])
        
        bin_counts, _, patches = (plt if ax == None else ax).hist(stat_values, bin_edges, 
                                                                    alpha=0.7,
                                                                    label=h_name,
                                                                    linestyle="solid",
                                                                    linewidth=1.5,
                                                                    hatch="///",
                                                                    color="w",
                                                                    histtype="step",
                                                                    ec=colorpalette[i],
                                                                    weights=np.ones_like(stat_values)/(len(stat_values) if normalize else 1))
        
        if isinstance(x_in, str):
            if callable(fit_func):
                fit_data = fit_func(bin_centers)
                
                if normalize:
                    fit_data = fit_data/fit_data.sum()
                    
                MSE = ((bin_counts - fit_data)**2).sum()*1/(len(bin_counts))
                RMSE = sqrt(MSE)
                
                ss_res = ((bin_counts - fit_data) ** 2).sum()
                ss_tot = (((bin_counts - np.mean(bin_counts)) ** 2)).sum()
                COE = 1 - (ss_res / ss_tot) # R^2
                
                text_rms = format_st(MSE) if not scientific_stats else f'{MSE:.2E}'
                text_rmse = format_st(RMSE) if not scientific_stats else f'{RMSE:.2E}'
                text_rsq = f'{COE:.2f}' if not scientific_stats else f'{COE:.2E}'
                
                (plt if ax == None else ax).plot(bin_centers, fit_data, color="red", alpha=0.7)
                fig.text(text_start_x, text_start_y - text_spacing_y*2*i,
                 f"Fit{fit_func.__name__ if not fit_func.__name__ == '<lambda>' else ''}\nMSE: {text_rms}\nRMSE: {text_rmse}\nR^2: {text_rsq}" , # + ("" if not isinstance(fit_opts, dict) else "\n".join("{0}:{1:.2f}".format(key, fit_opts[key]) for key in fit_opts.keys()))
                 #color=colorpalette[i],
                 bbox=dict(edgecolor="red", facecolor="w"),
                 fontsize='medium' if legendsize is None else legendsize,
                 horizontalalignment='right',
                 verticalalignment='top',
                 transform=ax.transAxes)
            
        mean = np.average(stat_values)
        std_dev = np.std(stat_values)
        
        mean_stat = f'{mean:.2E}' if scientific_stats else f'{mean:.2f}'
        std_dev_stat = f'{std_dev:.2E}' if scientific_stats else f'{std_dev:.2f}'
        
        fig.text(text_start_x, text_start_y - text_spacing_y*((2*i+1) if callable(fit_func) else i),
                f"{h_name}\nEntries: {len(values)}\nMean: {mean_stat}\nStd Dev: {std_dev_stat}",
                #color=colorpalette[i],
                bbox=dict(edgecolor=colorpalette[i], facecolor="w"),
                fontsize='medium' if legendsize is None else legendsize,
                horizontalalignment='right',
                verticalalignment='top',
                transform=ax.transAxes)
    
    #plt_obj.legend(loc='upper right', bbox_to_anchor=(1.1, 1.05))
    
    """
    if ax == None: 
        if title is not None:
            plt.title(title, fontdict = {'fontsize' : titlesize})
            
        if xlabel is not None:
            plt.xlabel(xlabel, fontsize=fontsize)
            
        if ylabel is not None:
            plt.ylabel(ylabel, fontsize=fontsize)
            
        plt.xscale(xscale)
        plt.yscale(yscale)
                    
        plt.xlim(xlim_view)
        plt.show()
        
        if ylim is not None:
            ax.ylim(ylim)
    else:
    """
    
    plot_styling(ax, ticksize_minor, ticksize_major, xscale, yscale, ylim, xlabel, ylabel, title, fontsize, titlesize)
    ax.set_xlim(left=xlim_view[0], right=xlim_view[1])
    
    return fig

def plot_styling(ax, ticksize_minor:int=10, ticksize_major:Optional[int]=None,
                 xscale:str="linear", yscale:str="linear",
                 ylim=None, xlabel:Optional[str] = None, ylabel:Optional[str]=None, title:Optional[str]=None,
                 fontsize:Optional[Union[str, int]]=14, titlesize:Union[int, str]=15):
    
    ax.tick_params(axis='both', which='major', labelsize=(ticksize_minor+2 if ticksize_major is None else ticksize_major))
    ax.tick_params(axis='both', which='minor', labelsize=ticksize_minor)
    
    if title is not None:
        ax.set_title(title, fontdict = {'fontsize' : titlesize})
        
    ax.set_xscale(xscale)
    ax.set_yscale(yscale)
    
    if xlabel is not None:
        ax.set_xlabel(xlabel, fontsize=fontsize)
        
    if ylabel is not None:
        ax.set_ylabel(ylabel, fontsize=fontsize)
    
    if ylim is not None:
        ax.set_ylim(ylim)
    
def plot_confusion(conf_mat, fontsize=12, ticksize_minor:int=10, ticksize_major:Optional[int]=None, title:Optional[str]=None, titlesize:Union[int, str]=15):
    import seaborn as sns
    
    if ticksize_major is None:
        ticksize_major = ticksize_minor + 2
    
    TP = conf_mat[0][0]
    TN = conf_mat[1][1]

    FN = conf_mat[0][1]
    FP = conf_mat[1][0]
    
    # Predicted
    PP = TP+FP
    PN = TN+FN
    
    # Actual
    P = TP + FN
    N = TN + FP
    
    annot = np.array([
        [f"{TP}\n(TPR: {TP/P*100:.2f}%)", f"{FN}\n(FNR: {FN/P*100:.2f}%)"],
        [f"{FP}\n(FPR: {FP/N*100:.2f}%)", f"{TN}\n(TNR: {TN/N*100:.2f}%)"],
    ])

    conf_mat = pd.DataFrame([
        [TP, FN],
        [FP, TN]
    ], index=["Sig", "Bkg"], columns=["Sig", "Bkg"])

    ax = sns.heatmap(conf_mat, annot=annot, fmt = '')
    
    #ax.set_xlabel("Predicted label", fontsize=fontsize)
    #ax.set_ylabel("Actual label", fontsize=fontsize)
    
    ax.tick_params(axis='both', which='major', labelsize=ticksize_major)
    ax.tick_params(axis='both', which='minor', labelsize=ticksize_minor)
    
    plot_styling(ax, ticksize_minor=ticksize_minor, ticksize_major=ticksize_major,
                 xscale="linear", yscale="linear", ylim=None,
                 xlabel="Predicted label", ylabel="Actual label",
                 title=title, fontsize=fontsize, titlesize=titlesize)
    
