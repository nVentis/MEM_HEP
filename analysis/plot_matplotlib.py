import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
import matplotlib.pylab as pylab
import pandas as pd
from matplotlib import rcParams as rcp
from typing import Optional, Union
from math import sqrt
#from scipy.stats import chisquare

rcp['hatch.linewidth'] = 0.5  # previous pdf hatch linewidth
rcp['font.family'] = 'monospace' # per default use monospace fonts

def format_st(f, t = 0.01):
    return f"{f:.2f}" if f >= t else f"<{t}"

def fontsize(fs):
    pylab.rcParams.update({
        'legend.fontsize': fs,
        'axes.labelsize': fs,
        'axes.titlesize': fs,
        'xtick.labelsize': fs,
        'ytick.labelsize': fs})

def plot_hist(data:Union[dict,pd.DataFrame], x:Optional[Union[str,list]], fit_func = None, fit_opts:Optional[dict] = None, labels=None, colorpalette=None, bins=128, xlim_binning=False, xlim:Optional[list] = None, ylim=None, xlabel:Optional[str] = None, ylabel:Optional[str]=None, units="", normalize=False, title:Optional[str] = "Likelihood-Analysis", ax=None, filter_nan:bool=False, text_start_x:float=0.965, text_start_y:float=0.97, text_spacing_y:float=0.23, xscale = "linear", yscale = "linear", fontsize:Optional[Union[str, int]]=14, legendsize = None, titlesize:Union[int, str]=15, ticksize_minor:int=10, ticksize_major=None):
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
        xlim_binning (bool, optional): _description_. Defaults to False.
        xlim (Optional[list], optional): _description_. Defaults to None.
        ylim (_type_, optional): _description_. Defaults to None.
        xlabel (Optional[str], optional): _description_. Defaults to None.
        ylabel (Optional[str], optional): _description_. Defaults to None.
        units (str, optional): _description_. Defaults to "".
        normalize (bool, optional): _description_. Defaults to False.
        title (Optional[str], optional): _description_. Defaults to "Likelihood-Analysis".
        ax (_type_, optional): _description_. Defaults to None.
        filter_nan (): Whether or not to filter out nan data; useful for DataFrames of unequal sizes. Defaults to False.
        text_start_x (float, optional): _description_. Defaults to 1.02.
        text_start_y (float, optional): _description_. Defaults to 1.05.
        text_spacing_y (float, optional): Height of each textbox . Defaults to 0.11.
        xscale (str, optional): _description_. Defaults to "linear".
        yscale (str, optional): _description_. Defaults to "linear".
        fontsize (Optional[str], optional): _description_. Defaults to None.
        ticksize_major (int, optional): ticksize_minor+2 if None. Defaults to None
    """
    
    if ax == None:
        fig, ax = plt.subplots()
        fig.set_dpi(100)
        fig.set_figwidth(8)
        
        fig.set_figheight(6)
    else:
        fig = plt.gcf()
    
    if colorpalette is None:
        #colorpalette = ["tab:blue", "tab:red", "y", "tab:pink", "tab:cyan", "tab:olive"]
        prop_cycle = plt.rcParams['axes.prop_cycle']
        colorpalette = prop_cycle.by_key()['color']
    
    # If data is one-dimensional, only plot this
    xlim_view = xlim
    
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
        if len(list(data.shape)) == 1:
            columns = [None] # In this case, data is assumed to contain just one column of data, which is to be histogrammed
            if xlim_view is None:
                xlim_view = [0.98*data.min(), 1.02*data.max()]
        else:
            columns = x
            if xlim_view is None:
                xlim_view = [0.98*data[x].min().min(), 1.02*data[x].max().max()]

    for i in range(len(columns)):
        column = columns[i]
        values = data if column is None else data[column]
        
        if filter_nan == True:
            values = values[~np.isnan(values)]
            
        # Limits
        min_val = xlim[0] if (xlim is not None and xlim_binning) else np.min(values)
        max_val = xlim[1] if (xlim is not None and xlim_binning) else np.max(values)
        
        bin_edges = np.linspace(min_val, max_val, bins + 1)
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
            
        # Additional behavior if only one column is to be plotted
        h_name  = (x if isinstance(x, str) else "Some Plot") if column is None else (column if labels is None else labels[i])
        
        bin_counts, _, patches = (plt if ax == None else ax).hist(values, bin_edges,
                 alpha=0.7,
                 label=h_name,
                 linestyle="solid",
                 linewidth=1.5,
                 hatch="///",
                 color="w",
                 histtype="step",
                 ec=colorpalette[i],
                 weights=np.ones_like(values)/(len(values) if normalize else 1))
        
        if isinstance(x, str):
            if callable(fit_func):
                fit_data = fit_func(bin_centers)
                
                if normalize:
                    fit_data = fit_data/fit_data.sum()
                    
                MSE = ((bin_counts - fit_data)**2).sum()*1/(len(bin_counts))
                RMSE = sqrt(MSE)
                
                ss_res = ((bin_counts - fit_data) ** 2).sum()
                ss_tot = (((bin_counts - np.mean(bin_counts)) ** 2)).sum()
                COE = 1 - (ss_res / ss_tot) # R^2
                
                (plt if ax == None else ax).plot(bin_centers, fit_data, color="red")
                fig.text(text_start_x, text_start_y - text_spacing_y*2*i,
                 f"Fit{fit_func.__name__ if not fit_func.__name__ == '<lambda>' else ''}\nMSE: {format_st(MSE)}\nRMSE: {format_st(RMSE)}\nR^2: {COE:.2f}" , # + ("" if not isinstance(fit_opts, dict) else "\n".join("{0}:{1:.2f}".format(key, fit_opts[key]) for key in fit_opts.keys()))
                 #color=colorpalette[i],
                 bbox=dict(edgecolor="red", facecolor="w"),
                 fontsize='medium' if legendsize is None else legendsize,
                 horizontalalignment='right',
                 verticalalignment='top',
                 transform=ax.transAxes)
        
        fig.text(text_start_x, text_start_y - text_spacing_y*((2*i+1) if callable(fit_func) else i),
                f"{h_name}\nEntries: {len(values)}\nMean: {np.average(values):.2f}\nStd Dev: {np.std(values):.2f}",
                #color=colorpalette[i],
                bbox=dict(edgecolor=colorpalette[i], facecolor="w"),
                fontsize='medium' if legendsize is None else legendsize,
                horizontalalignment='right',
                verticalalignment='top',
                transform=ax.transAxes)
    
    #plt_obj.legend(loc='upper right', bbox_to_anchor=(1.1, 1.05))
    
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
        ax.tick_params(axis='both', which='major', labelsize=(ticksize_minor+2 if ticksize_major is None else ticksize_major))
        ax.tick_params(axis='both', which='minor', labelsize=ticksize_minor)
        
        if title is not None:
            ax.set_title(title, fontdict = {'fontsize' : titlesize})
            
        ax.set_xlim(xlim_view)
        ax.set_xscale(xscale)
        ax.set_yscale(yscale)
        
        if xlabel is not None:
            ax.set_xlabel(xlabel, fontsize=fontsize)
            
        if ylabel is not None:
            ax.set_ylabel(ylabel, fontsize=fontsize)
        
        if ylim is not None:
            ax.set_ylim(ylim)


def plot_confusion(conf_mat, normalize=False, fontsize=12, ticksize_minor:int=10, ticksize_major=None):
    import seaborn as sns
    
    if ticksize_major is None:
        ticksize_major = ticksize_minor + 2
    
    TP = conf_mat[0][0]
    TN = conf_mat[1][1]

    FN = conf_mat[1][0]
    FP = conf_mat[0][1]
    
    tot = TP+TN+FN+FP
    
    annot = np.array([
        [f"{TP} ({TP/tot*100:.2f}%)", f"{FP} ({FP/tot*100:.2f}%)"],
        [f"{FN} ({FN/tot*100:.2f}%)", f"{TN} ({TN/tot*100:.2f}%)"],
    ])
    
    if normalize is True:
        annot = np.array([
            [f"{TP/tot*100:.2f}% ({TP})", f"{FP/tot*100:.2f}% ({FP})"],
            [f"{FN/tot*100:.2f}% ({FN})", f"{TN/tot*100:.2f}% ({TN})"],
        ])
         
        TP = TP/tot
        TN = TN/tot
        
        FN = FN/tot
        FP = FP/tot

    conf_mat = pd.DataFrame([
        [TP, FP],
        [FN, TN]
    ], index=["Sig", "Bkg"], columns=["Sig", "Bkg"])

    ax = sns.heatmap(conf_mat, annot=annot, fmt = '')
    
    ax.set_xlabel("True label", fontsize=fontsize)
    ax.set_ylabel("Predicted label", fontsize=fontsize)
    
    ax.tick_params(axis='both', which='major', labelsize=ticksize_major)
    ax.tick_params(axis='both', which='minor', labelsize=ticksize_minor)