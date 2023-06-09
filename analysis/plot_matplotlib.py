import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
from matplotlib import rcParams as rcp

rcp['hatch.linewidth'] = 0.5  # previous pdf hatch linewidth

def plot_hist(data, x, labels=None, colorpalette=None, bins=128, xlabel="", ylabel="", units="", normalize=False, title="Likelihood-Analysis"):
    
    # Autoscaling
    g_min = 0.98*data[x].min().min()
    g_max = 1.02*data[x].max().max()
    
    fig, ax = plt.subplots()
    fig.set_dpi(100)
    fig.set_figwidth(8)
    fig.set_figheight(6)
    
    if colorpalette is None:
        colorpalette = ["tab:red", "tab:blue", "y", "tab:pink", "tab:cyan", "tab:olive"]

    for i in range(0, len(x)):
        column = x[i]
        values = data[column]
        
        h_name  = labels[i] if labels is not None else column
        
        min_val = np.min(values)
        max_val = np.max(values)
        bin_centers = np.linspace(min_val, max_val, bins)
        
        plt.hist(values, bin_centers,
                 alpha=0.7,
                 label=h_name,
                 linestyle="solid",
                 linewidth=1.5,
                 hatch="///",
                 color="w",
                 histtype="step",
                 ec=colorpalette[i],
                 density=normalize == True)
        plt.text(1.05, 1. - 0.18*i,
                 h_name + "\nEntries: {0}\nMean: {1:.2f}\nStd Dev: {2:.2f}".format(len(values), np.average(values), np.std(values)),
                 color=colorpalette[i],
                 bbox=dict(edgecolor=colorpalette[i], facecolor="w"),
                 fontsize=10,
                 horizontalalignment='left',
                 verticalalignment='center',
                 transform=ax.transAxes)
    
    #plt.legend(loc='upper right', bbox_to_anchor=(1.1, 1.05))
    plt.title(title)
    plt.xlim(g_min, g_max)
    plt.show()