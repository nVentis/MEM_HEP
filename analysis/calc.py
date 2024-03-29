import numpy as np
import pandas as pd   
from typing import List

def calc_nll_llr_dtf_delta(df:pd.DataFrame, col_signal:str = "zhh_sigma", col_bkg:str = "zzh_sigma") -> pd.DataFrame:
    """Calculates signal/background negative log-likelihood and the log-likelihood-ratio assuming delta distributions as transfer functions

    Args:
        df (pd.DataFrame): _description_
        col_signal (str, optional): _description_. Defaults to "zhh_sigma".
        col_bkg (str, optional): _description_. Defaults to "zzh_sigma".

    Returns:
        pd.DataFrame: _description_
    """
     
    df["zhh_nll"] = -np.log(df[col_signal])
    df["zzh_nll"] = -np.log(df[col_bkg])
    df["llr"]     =  np.log(df[col_signal]/df[col_bkg])
    df.reset_index(drop=True, inplace=True)
    
    return df

def get_kinematics(data, true:bool, i:int, perm=[1,2,3,4]) -> List[float]:
    lep_key = "true_lep" if true else "lep"
    parton_key = "parton" if true else "jet"
    
    return [
        data[f"{lep_key}1_e"][i],
        data[f"{lep_key}1_px"][i],
        data[f"{lep_key}1_py"][i],
        data[f"{lep_key}1_pz"][i],
        
        data[f"{lep_key}2_e"][i],
        data[f"{lep_key}2_px"][i],
        data[f"{lep_key}2_py"][i],
        data[f"{lep_key}2_pz"][i],
        
        data[f"{parton_key}{perm[0]}_e"][i],
        data[f"{parton_key}{perm[0]}_px"][i],
        data[f"{parton_key}{perm[0]}_py"][i],
        data[f"{parton_key}{perm[0]}_pz"][i],
        
        data[f"{parton_key}{perm[1]}_e"][i],
        data[f"{parton_key}{perm[1]}_px"][i],
        data[f"{parton_key}{perm[1]}_py"][i],
        data[f"{parton_key}{perm[1]}_pz"][i],
        
        data[f"{parton_key}{perm[2]}_e"][i],
        data[f"{parton_key}{perm[2]}_px"][i],
        data[f"{parton_key}{perm[2]}_py"][i],
        data[f"{parton_key}{perm[2]}_pz"][i],
        
        data[f"{parton_key}{perm[3]}_e"][i],
        data[f"{parton_key}{perm[3]}_px"][i],
        data[f"{parton_key}{perm[3]}_py"][i],
        data[f"{parton_key}{perm[3]}_pz"][i],
    ]
    
def filter_for_tf(data, parton=True, jet=True, true_lep=True, lepton=True, parton1_pdg = 5, parton2_pdg = 5):
    if parton:
        data = data[((data["parton1_e"] > 0) & (data["parton2_e"] > 0) & (data["parton3_e"] > 0) & (data["parton4_e"] > 0) )]
    
    if jet:
        data = data[((data["jet1_e"] > 0) & (data["jet2_e"] > 0) & (data["jet3_e"] > 0) & (data["jet4_e"] > 0))]
    
    if true_lep:
        data = data[((data["true_lep1_e"] > 0) & (data["true_lep2_e"] > 0))]
        
    if lepton:
        data = data[((data["lep1_e"] > 0) & (data["lep2_e"] > 0))]
        
    if parton1_pdg is not None:
        data = data[(data["parton1_pdg"] == parton1_pdg)]
        
    if parton2_pdg is not None:
        data = data[(data["parton2_pdg"] == parton2_pdg)]
        
    data = data.copy()
    data.reset_index(drop=True, inplace=True)
    
    return data