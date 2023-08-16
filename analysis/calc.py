import numpy as np
import pandas as pd   

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

def normal_dist(x, mu, sigma):
    return 1/(sigma*np.sqrt(2*np.pi)) * np.exp(-1/2 * ((x-mu)/sigma)**2)

def dtf_dbgauss(coeffs, E_jet, E_part):
    result = 1
    
    for i in range(4): # len(coeffs)
        a, b = coeffs[i]
        E_j = E_jet[i]
        E_p = E_part[i]
        
        result *= dbgauss(a, b, E_jet, E_part)
        
    return result

def calc_nll_llr_dtf_dbgauss(df, coeff = None, col_signal = "zhh_sigma", col_bg = "zzh_sigma"):
    p4_true, p4_meas = extract_event_p4(df)
    
    def dtf(coeff, p4_meas, p4_int):
        ...
        
    def integrand(coeff, p4_meas, p4_int, me_func):
        return dtf(coeff, p4_meas, p4_int)*me_func(p4_int)
    
    # deltra distributions assumed for...
    # - lepton energies and
    # - jet solid angles  
    
    likelihood_sig = integrate(integrand, [])
    likelihood_bkg = integrate(integrand, [])
    
def get_kinematics(data, true:bool, i:int):
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
        
        data[f"{parton_key}1_e"][i],
        data[f"{parton_key}1_px"][i],
        data[f"{parton_key}1_py"][i],
        data[f"{parton_key}1_pz"][i],
        
        data[f"{parton_key}2_e"][i],
        data[f"{parton_key}2_px"][i],
        data[f"{parton_key}2_py"][i],
        data[f"{parton_key}2_pz"][i],
        
        data[f"{parton_key}3_e"][i],
        data[f"{parton_key}3_px"][i],
        data[f"{parton_key}3_py"][i],
        data[f"{parton_key}3_pz"][i],
        
        data[f"{parton_key}4_e"][i],
        data[f"{parton_key}4_px"][i],
        data[f"{parton_key}4_py"][i],
        data[f"{parton_key}4_pz"][i],
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