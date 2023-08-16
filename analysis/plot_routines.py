from typing import Union,Optional,Dict
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
import pandas as pd
from analysis.import_data import combine_columns, split_true_zhh_zzh
from analysis.plot_matplotlib import plot_hist

def plot_summary(data:pd.DataFrame, name: str, fig_path:Optional[str] = None):
    """Plots a summary about error occurences

    Args:
        data (_type_): _description_
        name (str): _description_
        fig_path (Optional[str], optional): _description_. Defaults to None.
    """
    
    import seaborn as sns
    
    fig, ((ax11, ax12)) = plt.subplots(1, 2, figsize=(9,6))
    fig.suptitle("Process summary: {}".format(name))

    fig11 = sns.countplot(data[data["is_zhh"] == 1], x="error_code", ax=ax11)
    fig11.set_title("{} / ZHH".format(name))

    fig12 = sns.countplot(data[data["is_zzh"] == 1], x="error_code", ax=ax12)
    fig12.set_title("{} / ZZH".format(name))
    
    if fig_path is not None:
        fig.savefig(fig_path)
    else:
        plt.show()
        
    plt.close(fig)

def plot_nll(data:pd.DataFrame, name:str, zhh_path:Optional[str] = None, zzh_path:Optional[str] = None, zhh_nll:str = "zhh_nll", zzh_nll:str = "zzh_nll"):
    """Plots the negative-log likelihood for true ZHH and true ZZH events

    Args:
        data (pd.DataFrame): ZHH or ZZH event data
        name (str): label
        zhh_path (Optional[str], optional): _description_. Defaults to None.
        zzh_path (Optional[str], optional): _description_. Defaults to None.
    """
    
    true_zhh, true_zzh = split_true_zhh_zzh(data)
    
    fig, ax = plt.subplots()
    plot_hist(true_zhh, x = [zhh_nll, zzh_nll], title="ZHH event data", normalize=True, labels=["ZHH {}".format(name), "ZZH".format(name)], xlabel="nll", ax=ax)
    
    if zhh_path is not None:
        fig.savefig(zhh_path)
    else:
        plt.show()
        
    plt.close(fig)
    
    fig, ax = plt.subplots()
    plot_hist(true_zzh, x = [zhh_nll, zzh_nll], title="ZZH event data", normalize=True, labels=["ZHH {}".format(name), "ZZH".format(name)], xlabel="nll", ax=ax)
    
    if zzh_path is not None:
        fig.savefig(zzh_path)
    else:
        plt.show()
        
    plt.close(fig)
    
    
def plot_llr(data:pd.DataFrame, name:str, fig_path:Optional[str] = None, llr_column:str="llr"):
    """Plots the Log-Likelihood-Ratio for signal/background (ZHH/ZZH)

    Args:
        data (pd.DataFrame): assumes DataFrame with existing llr column for both ZHH and ZZH events (defined as signal/background), by whatever means the likelihoods were calculated (see analysis.calc)
        name (str): label
        fig_path (Optional[str], optional): path where plot should be saved. if None, will attempt to plot in running python session. Defaults to None.
    """
    true_zhh, true_zzh = split_true_zhh_zzh(data)

    llr = {
        "zhh_llr": true_zhh[llr_column],
        "zzh_llr": true_zzh[llr_column]
    }
    
    fig, ax = plt.subplots()
    plot_hist(llr, x = ["zhh_llr", "zzh_llr"], labels=["ZHH event data", "ZZH event data"], title="LLR: {}".format(name), text_start_x= 0.26, normalize=True, xlabel="llr", ax=ax)
    
    if fig_path is not None:
        fig.savefig(fig_path)
    else:
        plt.show()
        
    plt.close(fig)
    

def plot_nll_llr_overview(data:Dict[str,pd.DataFrame], fig_path:Optional[str] = None):
    """Plots overview of me_nll and llr for the given dictionary of pd.DataFrames in one overview plot

    Args:
        data (Dict[str,pd.DataFrame]): dictionary with entries of structure name: calc_llr_dtf_delta(filter_data(mcp/reco/..._raw)), where raw was imported with import_data
        fig_path (Optional[str], optional): path where plot should be saved. if None, will attempt to plot in running python session. Defaults to None.
    """
    fig, axes = plt.subplots(len(data.keys()), 3, figsize=(6*len(data.keys()), 7*3))
    order = ["mcp", "tjt", "tjs", "tjmr", "reco"]
    for i in range(len(order)):
        name = order[i]
        df = data[name]
        
        ax1, ax2, ax3 = axes[i]
        
        true_zhh, true_zzh = split_true_zhh_zzh(df)
        
        plot_hist(true_zhh, x = ["zhh_nll", "zzh_nll"], title="ZHH event data", normalize=True, labels=["ZHH {}".format(name), "ZZH".format(name)], xlabel="nll", ax=ax1)
        plot_hist(true_zzh, x = ["zhh_nll", "zzh_nll"], title="ZZH event data", normalize=True, labels=["ZHH {}".format(name), "ZZH".format(name)], xlabel="nll", ax=ax2)
        
        llr = combine_columns({ "zhh_llr": true_zhh["llr"], "zzh_llr": true_zzh["llr"] })
        plot_hist(llr, x = ["zhh_llr", "zzh_llr"], labels=["ZHH event data", "ZZH event data"], title="LLR: {}".format(name), normalize=True, xlabel="llr", ax=ax3)
        
    if fig_path is not None:
        fig.savefig(fig_path)
    else:
        plt.show()
        
    plt.close(fig)