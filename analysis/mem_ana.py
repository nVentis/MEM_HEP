import pandas as pd
import numpy as np

from math import sqrt,pi,floor,log10
from analysis.mc.tools import variance_weighted_result
from typing import List, Optional, Union
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
    "system_pz": 0,
    "B_Z_bb": 0.1512, # 
    "B_H_bb": 0.569, # https://pdg.lbl.gov/2023/reviews/rpp2022-rev-higgs-boson.pdf  
    "pb_to_1oGeV2": 2.56819e-9
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

def parse_lines(lines, line_zhh, old:bool=False):
    result = [0., 0.]
    
    if old:
        if line_zhh > 0: result[0] = float(lines[line_zhh].split(": ")[1])
        if line_zhh > 0: result[1] = float(lines[line_zhh+1].split(": ")[1])
    else:
        if line_zhh > 0: result[0] = float(lines[line_zhh].split(": [")[1].split("]")[0])
        if line_zhh > 0: result[1] = float(lines[line_zhh+1].split(": [")[1].split("]")[0])
        
    return result

def get_result(event_dir:str, event_idx:int)->List[float]:
    #print(f"{event_dir}/event_{str(event_idx)}/result.txt")
    with open(f"{event_dir}/event_{str(event_idx)}/result.txt", "r") as f:
        lines = f.readlines()
        line_zhh = 0
        for line in lines:
            if line.startswith("ZHH:"):
                break
            
            line_zhh += 1
        
        result = parse_lines(lines, line_zhh, old=('[' not in lines[line_zhh]))
    
    return result

def get_result_npy(path:str,
                   perms_all:bool=True, perm_list:List[int]=[0])->Union[List[float],None]:

    f = np.load(path, allow_pickle=True).item()
    
    zhh_means = np.array(f['sig']['means'])
    zzh_means = np.array(f['bkg']['means'])
    
    if perms_all:
        perm_list = np.arange(len(zhh_means))
        
    zhh_means = zhh_means[perm_list]
    zzh_means = zzh_means[perm_list]
    
    if len(zhh_means) != len(zzh_means):
        return None
    
    return [zhh_means.mean(), zzh_means.mean()]

def sig_to_bkg(zhh_cross_sec:float=constants["sigma_zhh"], zzh_cross_sec:float=constants["sigma_zzh"],
               z_bb_branching:float=constants["B_Z_bb"], h_bb_branching:float=constants["B_H_bb"],)->float:
    
    return (zhh_cross_sec*h_bb_branching)/(zzh_cross_sec*z_bb_branching)

def finalize_mems(results:np.ndarray,
                  assume_zzh:bool=False, assume_zhh:bool=False,
                  pb_to_1oGeV2:float=constants["pb_to_1oGeV2"],
                  prefac:float=1/((2**24)*(pi**18)),
                  zhh_cross_sec:float=constants["sigma_zhh"], zzh_cross_sec:float=constants["sigma_zzh"],
                  z_bb_branching:float=constants["B_Z_bb"], h_bb_branching:float=constants["B_H_bb"]):
    
    if not assume_zhh and not assume_zzh: raise Exception('Invalid state')
        
    return results*prefac/(zhh_cross_sec if assume_zhh else zzh_cross_sec)*h_bb_branching*(h_bb_branching if assume_zhh else z_bb_branching)*pb_to_1oGeV2

def load_results(event_dir:str, reco:pd.DataFrame,
                 perms_all:bool=True, perm_list:List[int]=[0], use_npy:Optional[bool]=None,
                 normalize_samples:bool=True, pb_to_1oGeV2:float=constants["pb_to_1oGeV2"], add_generator:bool=False,) -> pd.DataFrame:
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
    import os.path as osp
    from glob import glob
    
    prefac = 1/((2**24)*(pi**18))
    
    events = np.array([int(name.replace("event_", "")) for name in listdir(event_dir)])
    results = []
    mask = np.zeros(len(events), dtype=bool)

    for i in range(len(events)):
        event_idx = events[i]
        if (use_npy != True) and perms_all:
            if osp.isfile(f"{event_dir}/event_{str(event_idx)}/result.txt"):
                mask[i] = True
                results.append(get_result(event_dir, event_idx))
        else:
            files = glob(f"{event_dir}/event_{str(event_idx)}/summary*.npy")
            if len(files):
                res = get_result_npy(path=files[0], perms_all=perms_all, perm_list=perm_list)
                if res is not None:
                    mask[i] = True
                    results.append(res)
                else:
                    print(event_idx)
            else:
                print(f'No file for event {event_idx}')
        
    results = np.array(results).T
    events = events[mask]
    
    assert(len(events) == results.shape[1])
    
    results = {
        "event_idx": events,
        "event": reco.iloc[events]['event'],
        "zhh_mem": finalize_mems(results[0], assume_zhh=True),
        "zzh_mem": finalize_mems(results[1], assume_zzh=True),
        "is_zhh": reco["is_zhh"][events],
        "is_zzh": reco["is_zzh"][events]
    }
    if add_generator:
        results["zhh_true"] = reco["zhh_mg5"][events]
        results["zzh_true"] = reco["zzh_mg5"][events]

    results = pd.DataFrame(results)
    results = results[((results["zhh_mem"] > 0) & (results["zzh_mem"] > 0))]
    
    if normalize_samples:
        stob = sig_to_bkg()
        btos = stob**-1
        
        sig_size = np.count_nonzero(results["is_zhh"])
        bkg_size = np.count_nonzero(results["is_zzh"])
        
        print('sig', sig_size, 'bkg', stob*bkg_size)
        
        if abs(sig_size - stob*bkg_size) > 2:
            if bkg_size < btos*sig_size:
                sig_size_max = round(stob*bkg_size)
                
                # Use complete background sample, lower signal fraction
                idx = np.where(results["is_zhh"] == 1)[0]
                np.random.seed(42)
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
    
def best_threshold(results, vals=None, r_column="r",
                   return_df=False, optimization_scheme:int=0, nsteps:int=1000):
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
        vals = np.linspace(np.min(results[r_column]), np.max(results[r_column]), nsteps)
    
    best = 9999
    best_t = np.max(results[r_column])
    
    stob = sig_to_bkg()
    
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
            cur = sqrt((TP - stob*TN)**2) + FN + FP #+ sqrt((FN-FP)**2)
        elif optimization_scheme == 1: # variation of 0
            cur = sqrt((TP - stob*TN)**2) + sqrt(FN**2 + FP**2)
        elif optimization_scheme == 2: # signal/background-ratio as expected, high TPR
            cur = sqrt(((TP - stob*TN)/P)**2) - TPR
        elif optimization_scheme == 3: # high TPR, high TNR 
            cur = -TPR -TNR 
        elif optimization_scheme == 4: # ratio of predicted positive/negative = stob
            if PN == 0:
                cur = 99999
            else:
                cur = (PP/PN - stob)**2            
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
        return best_t, pd.DataFrame(result)
    else:
        return best_t

def plot_r(zhh_data:Optional[np.ndarray]=None,
           zzh_data:Optional[np.ndarray]=None,
           name:Optional[str]=None, text_start_y:float=0.95, text_start_x:float=0.93, normalize=True, bins=64, kwargs:dict={}):
    from analysis.plot_matplotlib import plot_hist
    from analysis.import_data import split_true_zhh_zzh
    from matplotlib import pyplot as plt
    
    llr = {}
    labels = []
    
    if zhh_data is not None:
        llr['zhh_r'] = zhh_data
        labels.append('ZHH Events')
        
    if zzh_data is not None:
        llr['zzh_r'] = zzh_data
        labels.append('ZZH Events')

    args = { 'labels': labels, 'title': r"$D_{sig}$"+(f' {name}' if name is not None else ''),
                     'text_start_y': text_start_y, 'text_start_x': text_start_x, 'normalize': normalize,
                     'xlabel': r"$D_{sig}$", 'bins': bins, 'xscale': 'log', 'scientific_stats': True }
    
    plot_args = { **args, **kwargs }
    
    return plot_hist(llr, **plot_args)

