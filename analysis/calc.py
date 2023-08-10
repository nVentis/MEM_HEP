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

def extract_event_p4(df):
    return ([], [])
    
def me_sig(p4):
    ...
    
def me_bkg(p4):
    ...
    
def integrate(some_func, bounds):
    ...

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
    
    