import numpy as np

def calc_nll(df, col_signal = "zhh_sigma", col_bg = "zzh_sigma", col_signal_nll = "zhh_nll", col_bg_nll = "zzh_nll"):
    df[col_signal_nll] = -np.log(df[col_signal])
    df[col_bg_nll]     = -np.log(df[col_bg])
    df["llr"]          =  np.log(df[col_signal]/df[col_bg])
    df.reset_index(drop=True, inplace=True)
    
    return df