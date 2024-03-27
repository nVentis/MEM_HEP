import pandas as pd
import numpy as np

from math import sqrt,pi,floor,log10
from typing import List
import re

constants = {
    "m_b": 4.8, # truth MC: 4.8(?); wikipedia: 4.18
    "m_H": 125.,
    "m_Z": 91.19,
    "m_mu": 0.1056357046473643, # from average over all truth MC muon/anti-muons; wikipedia: 0.105658
    "sqrt_s": 500.,
    "sigma_zhh": 0.0111658, # from dumpevent /pnfs/desy.de/ilc/prod/ilc/mc-2020/ild/dst-merged/500-TDR_ws/hh/ILD_l5_o1_v02_nobg/v02-02-03/00015806/000/rv02-02-03.sv02-02-03.mILD_l5_o1_v02_nobg.E500-TDR_ws.I403001.Pe2e2hh.eL.pR.n000.d_dstm_15806_0.slcio
    "sigma_zzh": 0.0822856, # from dumpevent /pnfs/desy.de/ilc/prod/ilc/mc-2020/ild/dst-merged/500-TDR_ws/hh/ILD_l5_o1_v02_nobg/v02-02-03/00015806/000/rv02-02-03.sv02-02-03.mILD_l5_o1_v02_nobg.E500-TDR_ws.I403011.Pe2e2qqh.eL.pR.n000.d_dstm_15806_0.slcio
    "system_px": 0,
    "system_py": 0,
    "system_pz": 0
}

def parse_line_to_float(line:str, with_uncert:bool=False) -> float:   
    main = re.sub('\(\S+\)', '', line)    
    main = float(re.sub('\S+: ', '', main))
    
    #print(line)
    if with_uncert:
        uncert = float(re.findall(r'\((\d+)\)', line)[0])
        uncert = uncert*10**(floor(log10(main)) -2)
        
        return main, uncert
    else:
        return main

def get_result(event_dir:str, event_idx:int)->List[float]:
    #print(f"{event_dir}/event_{str(event_idx)}/result.txt")
    f = open(f"{event_dir}/event_{str(event_idx)}/result.txt", "r")
    lines = f.readlines()
    line_zhh = 0
    for line in lines:
        if line.startswith("ZHH:"):
            break
        
        line_zhh += 1
    
    result = [0., 0.]
    if line_zhh > 0:
        result[0] = float(lines[line_zhh].split(": ")[1])
        
    if line_zhh > 0:
        result[1] = float(lines[line_zhh+1].split(": ")[1])
    
    return result

def load_results(event_dir:str, reco:pd.DataFrame, zhh_cross_sec:float=constants["sigma_zhh"], zzh_cross_sec:float=constants["sigma_zzh"],
                 normalize_samples:bool=True, pb_to_1oGeV2:float = 2.56819e-9, add_generator:bool=False) -> pd.DataFrame:
    """_summary_

    Args:
        event_dir (str): _description_
        reco (pd.DataFrame): _description_
        zhh_cross_sec (float, optional): _description_. Defaults to constants["sigma_zhh"].
        zzh_cross_sec (float, optional): _description_. Defaults to constants["sigma_zzh"].
        normalize_samples (bool, optional): Makes sure the sample sizes for zhh and zzh are normalized by cross section. Defaults to True.
        pb_to_1oGeV2 (float, optional): _description_. Defaults to 2.56819e-9.

    Returns:
        pd.DataFrame: _description_
    """
    from os import listdir
    
    prefac = 1/((2**24)*(pi**18))
    
    events = np.array([int(name.replace("event_", "")) for name in listdir(event_dir)])
    results = np.array([get_result(event_dir, event_idx) for event_idx in events]).T
    
    assert(len(events) == results.shape[1])
    
    results = {
        "event_idx": events,
        "event": reco.iloc[events]['event'],
        "zhh_mem": results[0]*prefac/(zhh_cross_sec*pb_to_1oGeV2),
        "zzh_mem": results[1]*prefac/(zzh_cross_sec*pb_to_1oGeV2),
        "is_zhh": reco["is_zhh"][events],
        "is_zzh": reco["is_zzh"][events]
    }
    if add_generator:
        results["zhh_true"] = reco["zhh_mg5"][events]
        results["zzh_true"] = reco["zzh_mg5"][events]

    results = pd.DataFrame(results)
    if False:
        results = results[((results["zhh_mem"] > 0) & (results["zzh_mem"] > 0))]
        
        if normalize_samples:
            sig_to_bkg = zhh_cross_sec/zzh_cross_sec # ca. 0.1
            bkg_to_sig = sig_to_bkg**-1 # ca. 10
            
            sig_size = np.count_nonzero(results["is_zhh"])
            bkg_size = np.count_nonzero(results["is_zzh"])
            
            if abs(sig_size - sig_to_bkg*bkg_size) > 2:
                if bkg_size < bkg_to_sig*sig_size:
                    sig_size_max = round(sig_to_bkg*bkg_size)
                    
                    # Use complete background sample, lower signal fraction
                    idx = np.where(results["is_zhh"] == 1)[0]
                    mask = np.random.choice(range(0, len(idx)), size=(sig_size-sig_size_max), replace=False)

                    results.drop(results.index[idx[mask]], inplace=True)
                else:
                    # Use complete signal sample, lower background fraction
                    raise Exception("Not implemented")
                    

        results["r"] = results["zhh_mem"]/(results["zhh_mem"] + results["zzh_mem"])
    
    return results

def conf_mat(results, threshold:float):
    TP = np.count_nonzero((results["is_zhh"]) & (results["r"] > threshold))
    TN = np.count_nonzero((results["is_zzh"]) & (results["r"] < threshold))

    FN = np.count_nonzero((results["is_zhh"]) & (results["r"] < threshold))
    FP = np.count_nonzero((results["is_zzh"]) & (results["r"] > threshold))
    
    return [
        [TP, FN],
        [FP, TN]
    ]
    
def best_threshold(results, vals=None, r_column="r", zhh_cross_sec:float=constants["sigma_zhh"], zzh_cross_sec:float=constants["sigma_zzh"], return_df=False, optimization_scheme:int=0):
    """_summary_

    Args:
        results (_type_): _description_
        vals (_type_, optional): _description_. Defaults to None.
        r_column (str, optional): _description_. Defaults to "r".
        zhh_cross_sec (float, optional): _description_. Defaults to constants["sigma_zhh"].
        zzh_cross_sec (float, optional): _description_. Defaults to constants["sigma_zzh"].
        return_df (bool, optional): _description_. Defaults to False.
        optimization_scheme (int, optional): 0: low FN and FP, bkg and sig by ratio of cross-section.

    Returns:
        _type_: _description_
    """
    
    if vals is None:
        vals = np.linspace(np.min(results[r_column]), np.max(results[r_column]), 1000)
    
    best = 9999
    best_t = np.max(results[r_column])
    
    sig_to_bkg = zhh_cross_sec/zzh_cross_sec
    
    result = {
        "TP": [],
        "TN": [],
        "FN": [],
        "FP": []
    }
    
    P = np.count_nonzero((results["is_zhh"]))
    N = np.count_nonzero((results["is_zzh"]))
    
    for thresh in vals:
        TP = np.count_nonzero((results["is_zhh"]) & (results[r_column] > thresh))
        TN = np.count_nonzero((results["is_zzh"]) & (results[r_column] < thresh))

        FN = np.count_nonzero((results["is_zhh"]) & (results[r_column] < thresh))
        FP = np.count_nonzero((results["is_zzh"]) & (results[r_column] > thresh))
        
        PP = TP+FP
        PN = TN+FN
        
        TPR = TP/P
        TNR = TN/N
        
        cur = 0
        if optimization_scheme == 0: # signal/background-ratio as expected
            cur = sqrt((TP - sig_to_bkg*TN)**2) + FN + FP #+ sqrt((FN-FP)**2)
        elif optimization_scheme == 1: # variation of 0
            cur = sqrt((TP - sig_to_bkg*TN)**2) + sqrt(FN**2 + FP**2)
        elif optimization_scheme == 2: # signal/background-ratio as expected, high TPR
            cur = sqrt(((TP - sig_to_bkg*TN)/P)**2) - TPR
        elif optimization_scheme == 3: # high TPR, high TNR 
            cur = -TPR -TNR 
        elif optimization_scheme == 4: # ratio of predicted positive/negative = sig_to_bkg
            if PN == 0:
                cur = 99999
            else:
                cur = (PP/PN - sig_to_bkg)**2            
        else:
            raise Exception("Unknown scheme")
            
        if cur < best:
            best = cur
            best_t = thresh
            
        result["TP"].append(TP)
        result["TN"].append(TN)
        result["FP"].append(FP)
        result["FN"].append(FN)
    
    if return_df:
        return pd.DataFrame(result)
    else:
        return best_t

def plot_r(data, name, yscale="log", text_start_y:float=0.95, text_start_x:float=0.93, normalize=True, bins=64):
    from analysis.plot_matplotlib import plot_hist
    from analysis.import_data import split_true_zhh_zzh
    from matplotlib import pyplot as plt
    
    true_zhh, true_zzh = split_true_zhh_zzh(data)

    llr = {
        "zhh_r": true_zhh["r"],
        "zzh_r": true_zzh["r"]
    }

    fig, ax = plt.subplots()
    plot_hist(llr, x = ["zhh_r", "zzh_r"], labels=["ZHH event data", "ZZH event data"], title=r"$D_{bkg}$ " + f"({name})", text_start_y=text_start_y, text_start_x=text_start_x, normalize=normalize, xlim=(-0.02,1.02), xlim_binning=(0,1.), xlabel=r"$D_{bkg}$", ax=ax, bins=bins, yscale=yscale)

