import pandas as pd
import numpy as np
import os.path as osp
import os
from pathlib import Path
from analysis.split_event_tree import ttype_column
from analysis.mem_ana import constants
from analysis.convert import root_to_numpy
from tqdm.auto import tqdm
from typing import Optional

cache_dir = "/nfs/dust/ilc/user/bliewert/fullflow_v3/comparison/cache"

def import_data(path):
    # Load numpy file to pandas dataframe
    df = pd.DataFrame(np.load(path, allow_pickle=True))
    
    # For total ME, sum over helicity final states (for Z boson with spin s=1, 2s+1=3 possibilities, s_z = -1,0,1), average over helicity initial states (RL, LR)
    # Here, however, we just want to summarize the 
    df["zzh_sigmalr"] = 1/2 * ( df["zzh_sigmalrl"] + df["zzh_sigmalrr"] )
    
    # Add a column called true_type which is either "zhh" or "zzh", based on the columns "is_zhh" and "is_zzh"
    ttype_column(df)
    
    return df

def samples_set_ratio(
    df:pd.DataFrame,
    sig_to_bkg:float,
    is_sig_col:str
):
    """Assumes len(df) = sig_size*(1+sig_to_bkg**-1) = sig_size + bkg_size

    Args:
        df (pd.DataFrame): _description_
        sig_to_bkg (float): ratio of signal/background
        is_sig_col (str): column in df that's 1 for sig, 0 for bkg

    Raises:
        Exception: _description_

    Returns:
        _type_: _description_
    """
    out = df
    sig_size = np.count_nonzero(df[is_sig_col])
    
    bkg_size = len(df) - sig_size
    sig_size_max = int(sig_to_bkg*bkg_size)
    
    # Use complete background sample, lower signal fraction
    idx = np.where(df[is_sig_col] == 1)[0]
    mask = np.random.choice(range(0, len(idx)), size=(sig_size-sig_size_max), replace=False)

    out = df.drop(df.index[idx[mask]])
        
    return out

def normalize_samples(
    df:pd.DataFrame,
    constants:dict=constants,
    comparison=["zhh", "zzh"]) -> pd.DataFrame:
    
    sig_cross_sec = constants[f"sigma_{comparison[0]}"]
    bkg_cross_sec = constants[f"sigma_{comparison[1]}"]
    
    sig_to_bkg = sig_cross_sec/bkg_cross_sec # ca. 0.1
    bkg_to_sig = sig_to_bkg**-1 # ca. 10
    
    sig_size = np.count_nonzero(df[f"is_{comparison[0]}"])
    bkg_size = np.count_nonzero(df[f"is_{comparison[1]}"])
    
    out = df
    if abs(sig_size - sig_to_bkg*bkg_size) > 2:
        if bkg_size < bkg_to_sig*sig_size:
            out = samples_set_ratio(out, sig_to_bkg=sig_to_bkg, is_sig_col=f"is_{comparison[0]}")
        else:
            # Use complete signal sample, lower background fraction
            raise Exception("Not implemented")
        
    return out

def import_true_reco(
    comparison=["zhh", "zzh"],
    src_file:Optional[str] = None,
    src_dir:str = "/nfs/dust/ilc/user/bliewert/fullflow_v3/comparison/",
    tree_name:str="dataTree",
    file_name:str="compare_reco.root",
    normalize:bool=False,
    equal_size:bool=False,
    constants:dict=constants,
    b1_decay_pdg:Optional[int]=5,
    b2_decay_pdg:Optional[int]=5,
    use_cache:bool=True,
    recalc:bool=False) -> pd.DataFrame:
    
    """Combines the raw ROOT data to a pandas dataframe. Optionally, ensures normalization according to cross-section
    Uses caching, i.e. saves the full sample as numpy array/pandas DataFrame when it's accessed the first time 

    Args:
        comparison (list, optional): Two-element list with signal first and background second. Defaults to ["zhh", "zzh"].
        src_file (str, optional): if supplied, events df will be loaded from here. otherwise either from cache or comparison directory
        src_dir (str, optional): _description_. Defaults to "/nfs/dust/ilc/user/bliewert/fullflow_v3/comparison/".
        tree_name (str, optional): _description_. Defaults to "dataTree".
        file_name (str, optional): _description_. Defaults to "compare_reco.root".
        normalize (bool, optional): _description_. Defaults to True.
        equal_size (bool, optional): _description_. Defaults to False.
        constants (dict, optional): _description_. Defaults to constants. Must contain cross sections sigma_{comparison[0]} and sigma_{comparison[1]}
        b1_decay_pdg (Optional[int], optional): _description_. Defaults to 5 (bottom).
        b2_decay_pdg (Optional[int], optional): _description_. Defaults to 5 (bottom).

    Raises:
        Exception: _description_

    Returns:
        pd.DataFrame: _description_
    """
    
    from os import listdir
    import os.path as osp
    Path(f'{src_dir}/cache').mkdir(parents=True, exist_ok=True)
    
    results = listdir(src_dir)
    results.sort()
    results.remove("log")
    results.remove("cache")
    
    df = None
    serialized_name = f'comparison_reco_{comparison[0]}_{comparison[1]}.npy'
    cache_path = osp.join(cache_dir, serialized_name)
    
    if use_cache and recalc and osp.isfile(cache_path):
        os.remove(cache_path)
    
    if src_file is not None:
        df = np.load(src_file, allow_pickle=True)
    elif use_cache == False or not osp.isfile(cache_path):
        df = root_to_numpy(osp.join(src_dir, results[0], "root/prod", file_name), tree_name, null_on_not_found=True)
        for i in (pbar := tqdm(range(1, len(results)))):
            pbar.set_description(f'Current length: {len(df)}')
            part_path = osp.join(src_dir, results[i], "root/prod", file_name)
            if osp.isfile(part_path):
                res = root_to_numpy(part_path, tree_name, null_on_not_found=True)
                if res is not None:
                    df = np.concatenate((df, res))
        
        if use_cache:
            np.save(cache_path, df, allow_pickle=True)
            print(f'Saved cache file to {cache_path}')
    else:
        print(f'Using cached file from {cache_path}')
        df = np.load(cache_path, allow_pickle=True)
        
    df = pd.DataFrame(df)
    
    if b1_decay_pdg is not None:
        df = df[(df["true_h1_decay_pdg"] == b1_decay_pdg)]
        
    if b2_decay_pdg is not None:
        df = df[((df["true_h2_decay_pdg"] == b2_decay_pdg) | (df["true_z2_decay_pdg"] == b2_decay_pdg))]
       
    if equal_size:
        df = samples_set_ratio(df, 1, f"is_{comparison[0]}")
    elif normalize:
        df = normalize_samples(df, constants=constants, comparison=comparison)
    
    df.reset_index(drop=True, inplace=True)
    
    return df

def reco_and_target(df_in:pd.DataFrame, target_col:str = "is_zhh") -> pd.DataFrame:
    """Extracts a subset with reco kinematics and target for classification tasks

    Args:import_true_reco
        df_in (pd.DataFrame): _description_
        target_col (str, optional): _description_. Defaults to "is_zhh".

    Returns:
        _type_: _description_
    """
    feature_cols = []
    props = ["e", "px", "py", "pz"]
    for i in range(1, 5):
        for prop in props:
            feature_cols.append(f"jet{i}_{prop}")
            
            if i < 3:
                feature_cols.append(f"lep{i}_{prop}")

    features = df_in[feature_cols]
    target = df_in[target_col]

    task = pd.concat([features, target], axis=1)
    
    return task

# Filters out only events with ZHH/ZZH->µµbar+qqbar+qqbar AND without any error 
def filter_data(df, check_only_error=False):
    df = df[(df["error_code"] == 0)]
    if not check_only_error:
        df = df[((df["zhh_sigma"] > 0) & (df["zzh_sigma"] > 0))]
    
    df = df.copy()
    df.reset_index(drop=True, inplace=True)
    
    return df

def combine_columns(some_dict):    
    df_new = pd.DataFrame(some_dict).copy()
    df_new.reset_index(drop=True,inplace=True)
    
    return df_new

def split_true_zhh_zzh(df:pd.DataFrame):
    true_zzh = df[(df["is_zzh"] == 1)].copy()
    true_zhh = df[(df["is_zhh"] == 1)].copy()

    true_zzh.reset_index(drop=True, inplace=True)
    true_zhh.reset_index(drop=True, inplace=True)
    
    return (true_zhh, true_zzh)