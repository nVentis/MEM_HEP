import numpy as np
from math import pi, sqrt

fit_funcs = {
    "log_normal": lambda x, mu, sigma:
        1/(x*sigma*sqrt(2*pi))*np.exp(- ((np.log(x)-mu)**2)/(2*sigma**2)),
    "laplace": lambda x, mu, sigma:
        1/(2*sigma)*np.exp(-np.abs(x-mu)/sigma),
        
    "lorentz": lambda x, x0, gamma:
        1/(pi*gamma*(1 + ( (x-x0)/gamma )**2 )),
    
    "uniform": lambda x, a: a*np.ones(len(x)),
    
    "gauss": lambda x, x0, sigma:
        1/(sqrt(2*pi)*sigma)*np.exp(-(x-x0)**2/(2*sigma**2)),
    
    "dbgauss_old": lambda x, x0, sigma, a2, x02, sigma2:
        # x (x0, sigma, a2, x02, sigma2)
        1/(sqrt(2*pi)*(sigma + sigma2*a2**2))*(np.exp(-(x-x0)**2/(2*sigma**2)) + (a2**2)*np.exp(-(x-x02)**2/(2*sigma2**2)))
}

def fit_inits(df):
    return {
        "laplace": [np.average(df), 5],
        "lorentz": [np.average(df), 2],
        "uniform": [0.1],
        "gauss": [np.average(df), np.std(df)],
        "dbgauss_old": [0, np.std(df), 0.2, 0, 2*np.std(df)],
        "log_normal": [np.log(np.mean(df)) - np.sqrt(np.std(df))/2, np.std(df)]
    }