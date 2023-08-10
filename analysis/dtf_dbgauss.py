from scipy.optimize import curve_fit,minimize
from analysis.import_data import import_data,filter_data,combine_columns
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

def dbgauss(a, b, E_jet, E_part, JES):
    # Inspired by https://journals.aps.org/prd/abstract/10.1103/PhysRevD.74.092005
    
    """_summary_

    Args:
        a (_type_): _description_
        b (_type_): _description_
        E_jet (_type_): _description_
        E_part (_type_): _description_
        JES (float, optional): Jet energy scale. Defaults to 1.

    Returns:
        _type_: _description_
    """
    
    # Parametization according to thesis https://inspirehep.net/files/b8064e6fd21931e696b7b91410462128
    p1 = a[0] + b[0]*E_part
    p2 = a[1] + b[1]*E_part
    p3 = a[2] + b[2]*E_part
    p4 = a[3] + b[3]*E_part
    p5 = a[4] + b[4]*E_part
    # or: p = a + b*E_Part
    
    E_jet = E_jet/JES
    
    dE = E_jet - E_part
    
    return 1/(np.sqrt(2*np.pi)*(p2+p3*p5)) * (
             np.exp( -((dE-p1)**2)/(2*p2**2) ) +
        p3 * np.exp( -((dE-p4)**2)/(2*p5**2) )
    )
    
def dbgauss_wrapped(X, *args):
    E_jet, E_part = X
    a = [args[0], args[1], 0, args[2], args[3]]
    #a = args[0:5]
    b = args[4:9]
    JES = 1
    
    return dbgauss(a, b, E_jet, E_part, JES)    