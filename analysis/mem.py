import vegas, pickle
import numpy as np
from math import pi
from analysis.calc import get_kinematics

from typing import Dict
from math import sqrt,sin,cos,atan2,acos,pi
from itertools import product
from analysis.cffi.mg5.lib import lib_options

constants = {
    "m_b": 4.799957253499145, # truth MC: 4.8(?); wikipedia: 4.18
    "m_H": 125.,
    "m_Z": 91.19,
    "m_mu": 0.1056357046473643, # from average over all truth MC muon/anti-muons; wikipedia: 0.105658
    "sqrt_s": 500.,
    "sigma_sig_tot": 1., # calculated with MadGraph5 for ZHH
    "sigma_bkg_tot": 1., # calculated with MadGraph5 for ZZH
    "system_px": 0,
    "system_py": 0,
    "system_pz": 0
}

# Priors on Rhb and Thb
prior_args = (
    [77.17406236, 38.99397172],
    [1.57131127, 0.78054955]
)

# When using abs(), strictly speaking, its not a normal dist. anymore
Rhb_prior = lambda size: np.abs(np.random.normal(loc=prior_args[0][0], scale=prior_args[0][1], size=size))
Thb_prior = lambda size: np.random.normal(loc=prior_args[1][0], scale=prior_args[1][1], size=size)
Phb_prior = lambda size: np.random.uniform(-pi, pi, size)

def get_angles(x, y, z):
    theta = acos(z/sqrt(x**2 + y**2 + z**2))
    phi = atan2(y, x)
    
    return (theta, phi)

def unit_vec(theta, phi):
    return np.array([
        sin(theta)*cos(phi),
        sin(theta)*sin(phi),
        cos(theta)
    ])
    
def axis_vec(axis = 0) -> np.ndarray:
    return np.array([
        1. if axis == 0 else 0.,
        1. if axis == 1 else 0.,
        1. if axis == 2 else 0.
    ])
    
def get_kinematics_tf(data, idx:int, mode:int = 1, order = ["1", "2", "1", "2", "3", "4"]):
    """_summary_

    Args:
        data (_type_): event data
        idx (int): _description_
        mode (int, optional): 1 = Reco (default); 0 = MC Truth

    Returns:
        _type_: _description_
    """
    p = "jet" if mode == 1 else "parton"
    l = "lep" if mode == 1 else "true_lep"
    
    angles = np.array([
        get_angles(data[f"{l}{order[0]}_px"][idx], data[f"{l}{order[0]}_py"][idx], data[f"{l}{order[0]}_pz"][idx]),
        get_angles(data[f"{l}{order[1]}_px"][idx], data[f"{l}{order[1]}_py"][idx], data[f"{l}{order[1]}_pz"][idx]),
        get_angles(data[f"{p}{order[2]}_px"][idx], data[f"{p}{order[2]}_py"][idx], data[f"{p}{order[2]}_pz"][idx]),
        get_angles(data[f"{p}{order[3]}_px"][idx], data[f"{p}{order[3]}_py"][idx], data[f"{p}{order[3]}_pz"][idx]),
        get_angles(data[f"{p}{order[4]}_px"][idx], data[f"{p}{order[4]}_py"][idx], data[f"{p}{order[4]}_pz"][idx]),
        get_angles(data[f"{p}{order[5]}_px"][idx], data[f"{p}{order[5]}_py"][idx], data[f"{p}{order[5]}_pz"][idx])
    ])
    
    energies = np.array([
        data[f"{l}{order[0]}_e"][idx],
        data[f"{l}{order[1]}_e"][idx],
        data[f"{p}{order[2]}_e"][idx],
        data[f"{p}{order[3]}_e"][idx],
        data[f"{p}{order[4]}_e"][idx],
        data[f"{p}{order[5]}_e"][idx]
    ])
    
    momenta = np.array([
        sqrt(data[f"{l}{order[0]}_px"][idx]**2 + data[f"{l}{order[0]}_py"][idx]**2 + data[f"{l}{order[0]}_pz"][idx]**2),
        sqrt(data[f"{l}{order[1]}_px"][idx]**2 + data[f"{l}{order[1]}_py"][idx]**2 + data[f"{l}{order[1]}_pz"][idx]**2),
        sqrt(data[f"{p}{order[2]}_px"][idx]**2 + data[f"{p}{order[2]}_py"][idx]**2 + data[f"{p}{order[2]}_pz"][idx]**2),
        sqrt(data[f"{p}{order[3]}_px"][idx]**2 + data[f"{p}{order[3]}_py"][idx]**2 + data[f"{p}{order[3]}_pz"][idx]**2),
        sqrt(data[f"{p}{order[4]}_px"][idx]**2 + data[f"{p}{order[4]}_py"][idx]**2 + data[f"{p}{order[4]}_pz"][idx]**2),
        sqrt(data[f"{p}{order[5]}_px"][idx]**2 + data[f"{p}{order[5]}_py"][idx]**2 + data[f"{p}{order[5]}_pz"][idx]**2)        
    ])
    
    return (energies, momenta, angles)

def get_corrections(data, event_idx:int, is_reco:int = 1):
    qties = (
        [],
        [],
        [],
        []
    )
    for i, qty in enumerate(["e", "px", "py", "pz"]):
        for part in (["lep1", "lep2", "jet1", "jet2", "jet3", "jet4"] if is_reco else ["true_lep1", "true_lep2", "parton1", "parton2", "parton3", "parton4"]):
            qties[i].append(data[f"{part}_{qty}"][event_idx])
            
    return {
        "system_E": np.sum(qties[0]),
        "system_px": np.sum(qties[1]),
        "system_py": np.sum(qties[2]),
        "system_pz": np.sum(qties[3]),
    }

def get_parton_energies(data, event_idx):
    return [
        data["parton1_e"][event_idx],
        data["parton2_e"][event_idx],
        data["parton3_e"][event_idx],
        data["parton4_e"][event_idx]
    ]

def get_all_momenta(data, event_idx:int, as_4_vectors=False, is_reco=False):
    return np.concatenate(
        (get_lepton_momenta(data, event_idx, as_4_vectors, is_reco), get_parton_momenta(data, event_idx, as_4_vectors, is_reco))
    )
 
def get_parton_momenta(data, event_idx:int, as_4_vectors=False, is_reco=False):
    """_summary_

    Args:
        data (_type_): _description_
        event_idx (int): _description_
        as_4_vectors (bool, optional): If True, returns 4 vectors of partons, otherwise only momenta of 3 vectors. Defaults to False.

    Returns:
        _type_: _description_
    """
    
    p = "jet" if is_reco else "parton"
    
    if as_4_vectors:
        return [
            np.array([data[f"{p}1_e"][event_idx], data[f"{p}1_px"][event_idx], data[f"{p}1_py"][event_idx], data[f"{p}1_pz"][event_idx]]),
            np.array([data[f"{p}2_e"][event_idx], data[f"{p}2_px"][event_idx], data[f"{p}2_py"][event_idx], data[f"{p}2_pz"][event_idx]]),
            np.array([data[f"{p}3_e"][event_idx], data[f"{p}3_px"][event_idx], data[f"{p}3_py"][event_idx], data[f"{p}3_pz"][event_idx]]),
            np.array([data[f"{p}4_e"][event_idx], data[f"{p}4_px"][event_idx], data[f"{p}4_py"][event_idx], data[f"{p}4_pz"][event_idx]])
        ]
    else:
        return np.array([
            sqrt(data[f"{p}1_px"][event_idx]**2 + data[f"{p}1_py"][event_idx]**2 + data[f"{p}1_pz"][event_idx]**2),
            sqrt(data[f"{p}2_px"][event_idx]**2 + data[f"{p}2_py"][event_idx]**2 + data[f"{p}2_pz"][event_idx]**2),
            sqrt(data[f"{p}3_px"][event_idx]**2 + data[f"{p}3_py"][event_idx]**2 + data[f"{p}3_pz"][event_idx]**2),
            sqrt(data[f"{p}4_px"][event_idx]**2 + data[f"{p}4_py"][event_idx]**2 + data[f"{p}4_pz"][event_idx]**2)    
        ])
        
def get_lepton_momenta(data, event_idx:int, as_4_vectors=False, is_reco=False):
    """_summary_

    Args:
        data (_type_): _description_
        event_idx (int): _description_
        as_4_vectors (bool, optional): If True, returns 4 vectors of partons, otherwise only momenta of 3 vectors. Defaults to False.

    Returns:
        _type_: _description_
    """
    
    p = "lep" if is_reco else "true_lep"
    
    if as_4_vectors:
        return [
            np.array([data[f"{p}1_e"][event_idx], data[f"{p}1_px"][event_idx], data[f"{p}1_py"][event_idx], data[f"{p}1_pz"][event_idx]]),
            np.array([data[f"{p}2_e"][event_idx], data[f"{p}2_px"][event_idx], data[f"{p}2_py"][event_idx], data[f"{p}2_pz"][event_idx]]),
        ]
    else:
        return np.array([
            sqrt(data[f"{p}1_px"][event_idx]**2 + data[f"{p}1_py"][event_idx]**2 + data[f"{p}1_pz"][event_idx]**2),
            sqrt(data[f"{p}2_px"][event_idx]**2 + data[f"{p}2_py"][event_idx]**2 + data[f"{p}2_pz"][event_idx]**2),   
        ])
        
def prepare_int(data, event_idx:int, constants:dict, is_reco=True):
    """_summary_

    Args:
        data (_type_): _description_
        event_idx (int): _description_
        constants (dict): _description_
        is_reco (bool, optional): _description_. Defaults to True.

    Returns:
        _type_: _description_
    """
    
    p = "lep" if is_reco else "true_lep"
    
    mu1E = data[f"{p}1_e"][event_idx]
    mu1Th, mu1Ph = get_angles(data[f"{p}1_px"][event_idx], data[f"{p}1_py"][event_idx], data[f"{p}1_pz"][event_idx])
    mu1p = unit_vec(mu1Th, mu1Ph)*sqrt(mu1E**2 - constants["m_mu"]**2)
    
    mu2E = data[f"{p}2_e"][event_idx]
    mu2Th, mu2Ph = get_angles(data[f"{p}2_px"][event_idx], data[f"{p}2_py"][event_idx], data[f"{p}2_pz"][event_idx])
    mu2p = unit_vec(mu2Th, mu2Ph)*sqrt(mu2E**2 - constants["m_mu"]**2)
    
    return (
        (mu1E, mu1p),
        (mu2E, mu2p)
    )
    
def get_evt_constants(data, event_idx, constants=constants, use_reco=True):
    muon_kin = prepare_int(data, event_idx, constants, use_reco)
    ((mu1E, mu1p), (mu2E, mu2p)) = muon_kin
    get_parton_energies
    
    evt_constants = constants | {
        "system_E": constants["sqrt_s"] -mu1E -mu2E,
        "system_px": constants["system_px"] -mu1p[0] -mu2p[0],
        "system_py": constants["system_py"] -mu1p[1] -mu2p[1],
        "system_pz": constants["system_pz"] -mu1p[2] -mu2p[2],
    }
    
    return evt_constants

# INTEGRATION USING CFFI C++ BINDINGS

def int_bf_v2(data, event_idx:int, precond_size:int=4000000, mode:int=1, nitn:int=8, neval:int=16000000, nhcube_batch:int=100000, nwa:bool=lib_options["NWA"], use_tf:bool=True):
    """MEM integration in C++ using MG5 matrix elements

    Args:
        data (_type_): _description_
        event_idx (int): _description_
        precond_size (int, optional): _description_. Defaults to 4000000.
        mode (int, optional): _description_. Defaults to 1.
        nitn (int, optional): _description_. Defaults to 8.
        neval (int, optional): _description_. Defaults to 16000000.
        nhcube_batch (int, optional): _description_. Defaults to 100000.
        nwa (bool): whether or not to use NWA
        use_tf (bool, optional): whether or not to use the detector transfer function; False for debugging

    Returns:
        _type_: _description_
    """

    from analysis.cffi.mg5.lib import mc_batch, mc_batch_sigma
    
    reco_kin = get_kinematics(data, False, event_idx)
    
    # phi: from -pi to pi
    # theta: from 0 to pi
    
    bounds = [
        [0, pi],
        
        [-pi, pi],
        [0, 200],
        
        [0, pi],
        [-pi, pi],
        
        [0, 200],
        [0, pi]
    ]
    
    map = vegas.AdaptiveMap(bounds)
    x_arr = [
        Thb_prior(precond_size), # Thb1
        
        Phb_prior(precond_size), # Phb1
        Rhb_prior(precond_size), # Rhb1
        
        Thb_prior(precond_size), # Thb1b
        Phb_prior(precond_size), # Phb1b
        
        Rhb_prior(precond_size), # Rhb2
        Thb_prior(precond_size), # Thb2
    ]
    
    if not nwa:
        bounds.insert(0, [124.9**2, 125.1**2])
        x_arr.insert(0, np.ones(precond_size)*125**2) #mB1pow2
    
    x = np.stack(x_arr).T
    
    @vegas.batchintegrand
    def f(vars):
        print(f"PS points given:found [{len(vars)}:", end="")
        res = []
        if use_tf:
            res = mc_batch(reco_kin, vars.flatten().tolist(), mode=mode)
        else:
            res = mc_batch_sigma(vars.flatten().tolist(), mode=mode)
        
        found = np.count_nonzero(res)
        print(f"{found}] ({(100*found/len(vars)):6.2f} %)")
        return res
    
    map.adapt_to_samples(x, f(x), nitn=5)
    
    integ = vegas.Integrator(map, alpha=0.1, nhcube_batch=nhcube_batch)
    result = integ(f, nitn=nitn, neval=neval, save=f"{'zhh' if mode else 'zzh'}_{event_idx}_v3.pkl")
    
    print(result.summary())
    print('result = %s    Q = %.2f' % (result, result.Q))
    
    return result


def int_test_bf_v2(data, event_idx:int, samples_per_dim:int=4, evt_constants:dict=constants, use_reco:bool=True, mode:int=1, nwa:bool=lib_options["NWA"]):   
    """Tests the implementation using a uniform grid along the 8 axes of integration variables with samples_per_dim grid points per axis

    Args:
        data (_type_): _description_
        event_idx (int): _description_
        samples_per_dim (int, optional): _description_. Defaults to 4.
        evt_constants (dict, optional): _description_. Defaults to constants.
        use_reco (bool, optional): _description_. Defaults to True.
        mode (int, optional): 1 for signal (ZHH), 0 for background (ZZH)

    Returns:
        _type_: _description_
    """
    
    from analysis.cffi.mg5.lib import mc_batch
    
    reco_kin = get_kinematics(data, not use_reco,event_idx)
    
    boundaries = [
        (0., pi), # Thb1
        
        (-pi, pi), # Phb1
        (0., 200.), # Rhb1
        
        (0., pi), # Thb1b
        (-pi, pi), # Phb1b
        
        (0., 200.), # Rhb2
        (0., pi) # Thb2
    ]
    
    if not nwa:
        boundaries.insert(0, (124.8**2, 125.2**2)) # mB1_pow2
        
    dims = 7 if nwa else 8
    
    samples = []
    for boundary in boundaries:
        samples.append(np.linspace(boundary[0], boundary[1], samples_per_dim))
    
    int_variables = np.zeros(dims*samples_per_dim**dims) # 8 points per iteration; samples_per_dim**8 iterations
    for i, comb in enumerate(product(range(samples_per_dim), repeat=dims)):
        row = [samples[j][comb[j]] for j in range(dims)]
        int_variables[i*dims:(i+1)*dims] = row
            
    int_variables = int_variables.flatten().tolist()

    result = mc_batch(reco_kin, int_variables, mode=mode)
    
    return result

def get_true_int_args(data, event_idx:int, constants:dict, nwa=lib_options["NWA"]):
    energies, momenta, angles = get_kinematics_tf(data, event_idx, 0)
    
    Thb1,Phb1 = angles[2]
    Rhb1 = momenta[2]
    
    Thb1b, Phb1b = angles[3]
    
    Rhb2 = momenta[4]
    Thb2 = angles[4][0]
    
    int_variables = Thb1, Phb1, Rhb1, Thb1b, Phb1b, Rhb2, Thb2
    evt_constants = get_evt_constants(data, event_idx, constants, False)
    evt_constants["system_E"] = np.sum(energies[2:6])
    
    kin = get_kinematics(data, True, event_idx)
    
    return int_variables, evt_constants, kin