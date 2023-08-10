import click

import matplotlib.pyplot as plt
import logging
logger = logging.getLogger('mem_hep')

from analysis.import_data import import_data, filter_data
from analysis.plot_matplotlib import plot_hist
from math import sqrt, pow, pi
from typing import Optional
from os import path as osp
from os import makedirs, remove
import numpy as np

def plot_energy_transfer(data, name, plot_save_dir:Optional[str] = None, fit = "gauss", true_label = "parton", reco_label = "jet"):
    from scipy.optimize import curve_fit
    
    fig, axes = plt.subplots(1, 4, figsize=(6*len(data),8))
    fig.suptitle(name + r": $E_{" + reco_label + r"}-E_{" + true_label + r"}$", fontsize=18)
    
    for i in range(1,1+len(data)):
        n_bins = 128
        df = data[i-1]
        
        y, bins = np.histogram(df, bins=n_bins)
        x = (bins[:-1] + bins[1:]) / 2
        
        def bw(x, N, x0, G):
            return N*G/((2*pi)*((x-x0)**2 + ((G**2)/4)))
        
        def gauss(x, a, x0, sigma):
            return a*np.exp(-(x-x0)**2/(2*sigma**2))
        
        fit_func = bw if fit == "bw" else gauss
        
        init_bw = [len(df), np.average(df), 1]
        init_gauss = [len(df), np.average(df), np.std(df)]
        
        popt, pcov = curve_fit(fit_func, x, y, p0 = (init_bw if fit_func == bw else init_gauss))

        plot_hist(df, f"{reco_label.title()} {i}", fit_func=lambda x: fit_func(x, *popt), fit_opts=popt, bins=n_bins, xlim=(-100,100), ax=axes[i-1], xlabel=r"$ΔE$ [GeV]", title=f"{reco_label.title()} {i}", normalize=True, yscale="linear")
        #sns.histplot(data["jet{}_e".format(i)] - data["parton{}_e".format(i)], bins=128, ax=axes[i-1]).set_title("Jet {}".format(i))
        
    if plot_save_dir is not None:
        plot_path = osp.join(plot_save_dir, f"{name}.png")
        if osp.isfile(plot_path):
            remove(plot_path)
        
        fig.savefig(plot_path)
    else:
        plt.show()
        
    plt.close(fig)

def plot_jet_energy_transfer_from_df(data, name = "", plot_save_dir:Optional[str] = None, fit = "gauss"):
    df = []
    for i in range(1, 5):
        df.append(data["jet{}_e".format(i)] - data["parton{}_e".format(i)])
    
    plot_energy_transfer(df , name, plot_save_dir=plot_save_dir, fit=fit)

def plot_jet_energy_transfer_from_file(reco_npy_path:str, plots_path:str):
    import pandas as pd
    
    reco = filter_data(import_data(reco_npy_path), check_only_error=True)
    
    # Transfer function analysis for reco
    reco = reco[((reco["parton1_e"] > 0) & (reco["parton2_e"] > 0) & (reco["parton3_e"] > 0) & (reco["parton4_e"] > 0) )]
    reco = reco[((   reco["jet1_e"] > 0) & (   reco["jet2_e"] > 0) & (   reco["jet3_e"] > 0) & (   reco["jet4_e"] > 0) )]
    
    if not osp.isdir(plots_path):
        makedirs(plots_path)
    
    plot_jet_energy_transfer_from_df(reco, "ZHH+ZZH", plots_path)
    plot_jet_energy_transfer_from_df(reco[reco["is_zhh"] == 1], "ZHH", plots_path)
    plot_jet_energy_transfer_from_df(reco[reco["is_zzh"] == 1], "ZZH", plots_path)

@click.command("jet_et_plot")
@click.argument("npy_path")
@click.argument("plot_path")
def plot_jet_energy_transfer_from_file_command(npy_path:str, plot_path:str):
    plot_jet_energy_transfer_from_file(npy_path, plot_path)