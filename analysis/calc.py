import numpy as np

def calc_nll(df, zhh_me = "zhh_sigma", zzh_me = "zzh_sigma"):
    df["zhh_nll"] = -np.log(df[zhh_me])
    df["zzh_nll"] = -np.log(df[zzh_me])
    df["llr"]     = np.log(df[zzh_me]/df[zhh_me])
    df.reset_index(drop=True, inplace=True)
    
    return df