import numpy as np
import pandas as pd

# See https://stackoverflow.com/questions/40554179/how-to-keep-column-names-when-converting-from-pandas-to-numpy
def pd_to_np(df:pd.DataFrame):
    dtyp = np.dtype(list(zip(df.dtypes.index, df.dtypes)))
    arr_ip = [tuple(i) for i in df.to_numpy()]
    
    return np.array(arr_ip, dtype=dtyp)
    