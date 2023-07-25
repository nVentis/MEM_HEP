import pandas as pd
import numpy as np
from analysis.split_event_tree import ttype_column

def import_data(path):
    # Load numpy file to pandas dataframe
    df = pd.DataFrame(np.load(path, allow_pickle=True))
    
    # For total ME, sum over helicity final states (for Z boson with spin s=1, 2s+1=3 possibilities, s_z = -1,0,1), average over helicity initial states (RL, LR)
    # Here, however, we just want to summarize the 
    df["zzh_sigmalr"] = 1/2 * ( df["zzh_sigmalrl"] + df["zzh_sigmalrr"] )
    
    # Add a column called true_type which is either "zhh" or "zzh", based on the columns "is_zhh" and "is_zzh"
    ttype_column(df)
    
    return df

# Filters out only events with ZHH/ZZH->µµbar+qqbar+qqbar AND without any error 
def filter_data(df):
    df = df[(df["error_code"] == 0) & (df["zhh_sigma"] > 0) & (df["zzh_sigma"] > 0)]
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