import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
from matplotlib import rcParams as rcp
from typing import Optional, Union

rcp['hatch.linewidth'] = 0.5  # previous pdf hatch linewidth
rcp['font.family'] = 'monospace' # per default use monospace fonts

def plot_hist(data, x:Union[str,list], labels=None, colorpalette=None, bins=128, xlabel="", ylabel="", units="", normalize=False, title="Likelihood-Analysis", ax=None, text_start_x=1.05, text_start_y=1.05):
    """If x is a 
    """
    
    if ax == None:
        fig, ax = plt.subplots()
        fig.set_dpi(100)
        fig.set_figwidth(8)
        
        fig.set_figheight(6)
    else:
        fig = plt.gcf()
    
    if colorpalette is None:
        colorpalette = ["tab:red", "tab:blue", "y", "tab:pink", "tab:cyan", "tab:olive"]
    
    # If data is one-dimensional, only plot this
    if len(list(data.shape)) == 1:
        columns = [None] # In this case, data is assumed to contain just one column of data, which is to be histogrammed
        g_min = 0.98*data.min()
        g_max = 1.02*data.max()
    else:
        columns = x
        g_min = 0.98*data[x].min().min()
        g_max = 1.02*data[x].max().max()

    for i in range(len(columns)):
        column = columns[i]
        values = data if column is None else data[column]
        
        h_name  = (x if isinstance(x, str) else "Some Plot") if column is None else (column if labels is None else labels[i])
        
        min_val = np.min(values)
        max_val = np.max(values)
        bin_centers = np.linspace(min_val, max_val, bins)
        
        (plt if ax == None else ax).hist(values, bin_centers,
                 alpha=0.7,
                 label=h_name,
                 linestyle="solid",
                 linewidth=1.5,
                 hatch="///",
                 color="w",
                 histtype="step",
                 ec=colorpalette[i],
                 density=normalize == True)
        fig.text(text_start_x, text_start_y - 0.18*i,
                 h_name + "\nEntries: {0}\nMean: {1:.2f}\nStd Dev: {2:.2f}".format(len(values), np.average(values), np.std(values)),
                 color=colorpalette[i],
                 bbox=dict(edgecolor=colorpalette[i], facecolor="w"),
                 fontsize=9,
                 horizontalalignment='right',
                 verticalalignment='center',
                 transform=ax.transAxes)
    
    #plt_obj.legend(loc='upper right', bbox_to_anchor=(1.1, 1.05))
    if ax == None: 
        plt.title(title)
        plt.xlim(g_min, g_max)
        plt.show()
    else:
        ax.set_title(title)
        ax.set_xlim(g_min, g_max)