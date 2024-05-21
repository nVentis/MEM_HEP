import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
import matplotlib.pylab as pylab
import pandas as pd
from analysis.calc import calc_FWHM
from matplotlib import rcParams as rcp
from matplotlib.backends.backend_pdf import PdfPages
from typing import Optional, Union, Callable, Dict, List
import inspect
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

def export_figures(filename, figs=None, dpi=200):
    pp = PdfPages(filename)
    if figs is None:
        figs = [plt.figure(n) for n in plt.get_fignums()]
        #print(len(figs))
        
    for fig in figs:
        fig.savefig(pp, format='pdf')
    pp.close()

def plot_hist(data:Union[dict,pd.DataFrame], x:Optional[Union[str,list]]=None,
              fit_func:Optional[Callable]=None, fit_opts:Optional[Union[List[dict], dict]]=None,
              labels:Optional[List[str]]=None, colorpalette=None, bins:int=128,
              same_bins:bool=True, xlim_binning:Optional[Union[list,tuple]]=None, xlim:Optional[Union[list,tuple]]=None, ylim=None,
              xlabel:Optional[str] = None, ylabel:Optional[str]=None,
              title:Optional[str]=None, ax=None,
              normalize=False, filter_nan:bool=False, unitx:Optional[str]=None,
              with_fwhm:bool=False,
              text_start_x:float=0.965, text_start_y:float=0.97, text_spacing_y:float=0.22,
              xscale:Literal['linear', 'log']='linear', yscale:Literal['linear', 'log']="linear",
              fontsize:Optional[Union[str, int]]=14, legendsize=None, titlesize:Union[int, str]=15,
              ticksize_minor:int=10, ticksize_major=None,
              figsize:tuple=(8,6), figdpi:int=100, scientific_stats:bool=False,
              force_df:bool=True):
    """_summary_
    
    text_spacing_y: 0.11 for high-res
    

    Args:
        data (Union[dict,pd.DataFrame]): Can be a pd.DataFrame or simple dict (preferrable for columns of unequal size)
        x (Optional[Union[str,list]]): one column or a list of columns in data to be histogrammed.
        fit_func (Optional[function], optional): only supported if one column is to be plotted, i.e. x is a string
        fit_opts (Optional[Union[List[dict], dict]], opional): only used for printing
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
        if force_df:
            data = pd.DataFrame(data)
        else:
            for key in data:
                data[key] = np.array(data[key])
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
        
    if not isinstance(data, dict) and len(list(data.shape)) == 1:
        columns = [None] # In this case, data is assumed to contain just one column of data, which is to be histogrammed
        if xlim_view is None:
            xlim_view = [0.98*data.min(), 1.02*data.max()]
    else:
        columns = [x] if isinstance(x, str) else x
        if xlim_view is None:
            if isinstance(data, dict):
                xlim_view = (np.min(data[list(data.keys())[0]]), np.max(data[list(data.keys())[0]]))
                for key, dt in data.items():
                    xlim_view = (
                        min(xlim_view[0], dt.min()),
                        max(xlim_view[0], dt.max()),
                    )
            else:
                xlim_view = [0.98*data[x].min().min(), 1.02*data[x].max().max()]
                
    # If same_bins=True, infer limits and impose xlim_binning
    if xlim_binning is None:
        if same_bins or isinstance(data, pd.DataFrame):
            xlim_binning = xlim_view
        else: # why this?
            xlim_binning = [np.min(np.min(data)), np.max(np.max(data))]
    
    # TODO:    
    nextbox_y = text_start_y
    
    for i in range(len(columns)):
        column = columns[i]
        values = data if column is None else data[column]
        fit_fn = None if fit_func is None else (fit_func if callable(fit_func) else (fit_func[i] if (isinstance(fit_func, list) and callable(fit_func[i])) else None))
        
        if filter_nan == True:
            mask = np.isnan(values)
            if np.sum(mask):
                print(f'Warning: {np.sum(mask)} values NaN')
            
            values = values[~mask]
            
        # Limits
        min_val = xlim_binning[0] if (xlim_binning is not None) else np.min(values)
        max_val = xlim_binning[1] if (xlim_binning is not None) else np.max(values)
        
        if xscale == 'log':
            log_start, log_stop = np.log(min_val), np.log(max_val)
            bin_edges = np.exp(np.linspace(log_start, log_stop, bins+2))
        else:
            if xscale != 'linear':
                print(f'Unsupported xscale {xscale}. Reverting to linear')
                
            bin_edges = np.linspace(min_val, max_val, num=bins+1)
            
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2           

        # Filter for x_lim_binning
        stat_values = values
        if xlim_binning is not None:
            stat_values = stat_values[(stat_values >= xlim_binning[0]) & (stat_values <= xlim_binning[1])]
        
        # Additional behavior if only one column is to be plotted
        h_name  = (x_in if isinstance(x_in, str) else "Data") if column is None else (column if labels is None else labels[i])
        
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
        
        
        
        if callable(fit_fn):
            # Allow different parameters per fit/data
            fit_opt = None if fit_opts is None else (fit_opts[i] if len(np.shape(fit_opts)) == 2 else fit_opts)
            
            if fit_opt is not None:   
                fit_data = fit_fn(bin_centers, *fit_opt)
            else:
                fit_data = fit_fn(bin_centers)
            
            if normalize:
                fit_data = fit_data/fit_data.sum()
            
            MSE = ((bin_counts - fit_data)**2).sum()*1/(len(bin_counts))
            RMSE = sqrt(MSE)
            
            ss_reg = ((bin_counts - fit_data) ** 2).sum()
            ss_tot = (((bin_counts - np.mean(bin_counts)) ** 2)).sum()
            Rsquared = 1 - (ss_reg / ss_tot) # R^2
            
            text_rms = format_st(MSE) if not scientific_stats else f'{MSE:.2E}'
            text_rmse = format_st(RMSE) if not scientific_stats else f'{RMSE:.2E}'
            text_rsq = f'{Rsquared:.2f}' if not scientific_stats else f'{Rsquared:.2E}'
            
            fit_name = f' {fit_fn.__name__}' if not fit_fn.__name__ == '<lambda>' else ''
            fit_text = f"Fit{fit_name}\nMSE: {text_rms}\nRMSE: {text_rmse}\nR^2: {text_rsq}" # + ("" if not isinstance(fit_opts, dict) else "\n".join("{0}:{1:.2f}".format(key, fit_opts[key]) for key in fit_opts.keys()))
            fit_args = inspect.getfullargspec(fit_fn).args
            
            if fit_opt is not None and len(fit_args) > 1:
                # First argument is x and skipped
                keys = inspect.getfullargspec(fit_fn).args
                for j, val in enumerate(fit_opt):
                    key = keys[j+1]
                    fit_text += f'\n{key} = ' + (f'{val:.2E}' if scientific_stats else f'{val:.2f}')
            
            (plt if ax == None else ax).plot(bin_centers, fit_data, color="red", alpha=0.7)
            
            fig.text(text_start_x, text_start_y - text_spacing_y*2*i,
                fit_text,
                #color=colorpalette[i],
                bbox=dict(edgecolor="red", facecolor="w"),
                fontsize=(pylab.rcParams['legend.fontsize'] if legendsize is None else legendsize),
                horizontalalignment='right',
                verticalalignment='top',
                transform=ax.transAxes)
            
        mean = np.average(stat_values)
        std_dev = np.std(stat_values)
        fwhm = None
        
        # Extra statistics
        extra_text = ''
        if with_fwhm:
            fwhm = calc_FWHM(bin_centers, bin_counts)
            extra_text += f'\nFWHM: ' + (f'{fwhm:.2E}' if scientific_stats else f'{fwhm:.2f}')
        
        mean_stat = (f'{mean:.2E}' if scientific_stats else f'{mean:.2f}')+(f' {unitx}' if unitx is not None else '')
        std_dev_stat = (f'{std_dev:.2E}' if scientific_stats else f'{std_dev:.2f}')+(f' {unitx}' if unitx is not None else '')
        
        fig.text(text_start_x, text_start_y - text_spacing_y*((2*i+1) if callable(fit_fn) else i),
                f"{h_name}\nEntries: {len(stat_values)}\nMean: {mean_stat}\nStd Dev: {std_dev_stat}{extra_text}",
                #color=colorpalette[i],
                bbox=dict(edgecolor=colorpalette[i], facecolor="w"),
                fontsize=(pylab.rcParams['legend.fontsize'] if legendsize is None else legendsize),
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

def defaultsize(ls, ts:Optional[int]=None, axs:Optional[int]=None):
    if ts is None: ts = ls
    if axs is None: axs = ls
    
    pylab.rcParams.update({
        'legend.fontsize': ls,
        'axes.labelsize': axs,
        'axes.titlesize': axs,
        'xtick.labelsize': ts,
        'ytick.labelsize': ts})


defaultsize(12)

def plot_styling(ax,
                 ticksize_minor:Optional[int]=None, ticksize_major:Optional[int]=None,
                 xscale:Optional[str]="linear", yscale:Optional[str]="linear",
                 ylim=None, xlabel:Optional[str] = None, ylabel:Optional[str]=None, title:Optional[str]=None,
                 fontsize:Union[str, int, None]=None,
                 titlesize:Union[int, str, None]=None,
                 ticks_left:Union[bool,List,None]=True, ticks_bottom:Union[bool,List,None]=True):
    
    if fontsize is None: fontsize = pylab.rcParams['axes.labelsize'] 
    if titlesize is None: titlesize = fontsize + 2
    
    ticksize_minor_y = pylab.rcParams['ytick.labelsize'] - 2 if ticksize_minor is None else ticksize_minor
    ticksize_major_y = pylab.rcParams['ytick.labelsize']     if ticksize_major is None else ticksize_major
    
    ticksize_minor_x = pylab.rcParams['xtick.labelsize'] - 2 if ticksize_minor is None else ticksize_minor
    ticksize_major_x = pylab.rcParams['xtick.labelsize']     if ticksize_major is None else ticksize_major
    
    if isinstance(ticks_left, list):
        ax.set_yticks(range(len(ticks_left)), ticks_left, fontsize=pylab.rcParams['ytick.labelsize'])
    elif ticks_left is not None:
        if ticks_left:
            ax.tick_params(axis='y', which='major', labelsize=ticksize_major_y)
            ax.tick_params(axis='y', which='minor', labelsize=ticksize_minor_y)
        else:
            ax.tick_params(left=False, right=False, labelleft=False, labelright=False)
    
    if isinstance(ticks_bottom, list):
        ax.set_xticks(range(len(ticks_bottom)), ticks_bottom, fontsize=pylab.rcParams['xtick.labelsize'])
    elif ticks_bottom is not None:
        if ticks_bottom:
            ax.tick_params(axis='x', which='major', labelsize=ticksize_major_x)
            ax.tick_params(axis='x', which='minor', labelsize=ticksize_minor_x)
        else:
            ax.tick_params(bottom=False, top=False, labelbottom=False, labeltop=False)
    
    if title is not None: ax.set_title(title, fontdict = {'fontsize' : titlesize})
    
    if xscale is not None: ax.set_xscale(xscale)
    if yscale is not None: ax.set_yscale(yscale)
    
    if xlabel is not None: ax.set_xlabel(xlabel, fontsize=fontsize)
    if ylabel is not None: ax.set_ylabel(ylabel, fontsize=fontsize)
    
    if ylim is not None: ax.set_ylim(ylim)

def plot_lego(x:Optional[np.ndarray], y:Optional[np.ndarray],
              hist:Optional[np.ndarray], xedges:Optional[np.ndarray], yedges:Optional[np.ndarray],
              nbins:int=32, nbinsx:Optional[int]=None, nbinsy:Optional[int]=None,
              xlabel:Optional[str] = None, ylabel:Optional[str]=None,
              title:Optional[str]=None, ax=None, zscale:str='linear'):
    
    if hist is None:
        if nbinsx is None: nbinsx = nbins
        if nbinsy is None: nbinsy = nbins
        
        fig = plt.figure()          #create a canvas, tell matplotlib it's 3d
        if ax is None:
            ax = fig.add_subplot(111, projection='3d')

        #make histogram stuff - set bins - I choose 20x20 because I have a lot of data
        hist, xedges, yedges = np.histogram2d(x, y, bins=(20,20))
        
    xpos, ypos = np.meshgrid(xedges[:-1]+xedges[1:], yedges[:-1]+yedges[1:])

    xpos = xpos.flatten()/2.
    ypos = ypos.flatten()/2.
    zpos = np.zeros_like (xpos)

    if zscale == 'log':
        hist[hist != 0] = np.log(hist[hist != 0])

    dx = xedges [1] - xedges [0]
    dy = yedges [1] - yedges [0]
    
    dz = hist.flatten()

    cmap = cm.get_cmap('jet') # Get desired colormap - you can change this!
    max_height = np.max(dz)   # get range of colorbars so we can normalize
    min_height = np.min(dz)
    # scale each z to [0,1], and get their rgb values
    rgba = [cmap((k-min_height)/max_height) for k in dz] 

    ax.bar3d(xpos, ypos, zpos, dx, dy, dz, color=rgba, zsort='average')
    if zscale == 'log':
        ax.set_zscale('symlog')
    
    if title is not None: plt.title(title)
    if xlabel is not None: plt.xlabel(xlabel)
    if ylabel is not None: plt.ylabel(ylabel)
    
    return fig

def plot_confusion(conf_mat, fontsize=12, ticksize_minor:int=10, ticksize_major:Optional[int]=None,
                   title:Optional[str]=None, titlesize:Union[int, str]=15,
                   labels=['ZHH', 'ZZH'],):
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
    ], index=[labels[0], labels[1]], columns=[labels[0], labels[1]])

    ax = sns.heatmap(conf_mat, annot=annot, fmt = '')
    
    plot_styling(ax, ticksize_minor=ticksize_minor, ticksize_major=ticksize_major,
                 xscale=None, yscale=None, ylim=None,
                 xlabel="Predicted label", ylabel="Actual label",
                 title=title, fontsize=fontsize, titlesize=titlesize,
                 ticks_left=None, ticks_bottom=None,)
    
    return ax.get_figure()
    
def plot_roc(thresh_df:pd.DataFrame):
    from sklearn.metrics import auc
    
    TP = thresh_df["TP"]
    TN = thresh_df["TN"]
    
    FN = thresh_df["FN"]
    FP = thresh_df["FP"]
    
    P = TP + FN
    N = TN + FP
    
    FPR = FP/N
    TPR = TP/P

    fig, ax = plt.subplots()
    ax.plot(FPR, TPR)
    
    AUC = auc(FPR, TPR) #np.sum(TPR)/len(TPR)
    
    plot_styling(ax, title=f"ROC-curve [AUC:{AUC:.3f}]", xlabel="FPR", ylabel="TPR")
    
    return fig