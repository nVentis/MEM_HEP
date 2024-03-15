import vegas, pickle
import numpy as np
from math import sqrt, sin, cos, pi
from analysis.fit_funcs import fit_funcs
from analysis.mem import prepare_int, get_evt_constants, unit_vec, get_kinematics_tf, Phb_prior, Thb_prior, Rhb_prior

# Python (reference) implementation for phase space integration
# Does everything in Python, except evaluation of the matrix element

tf_E = lambda x: fit_funcs["lorentz"](x, *[-1.12594285,  6.24937028])
tf_Th = lambda x: fit_funcs["lorentz"](x, *[-3.90908961e-05,  1.96662831e-02])
tf_Ph = lambda x: fit_funcs["lorentz"](x, *[0.0001748,  0.02819419])

from math import acos,asin,copysign,atan2,atan,pi
from itertools import product

err_map = [-1, -2, -3, -4, -5, -6, -7, -8, -9, -10, -11]

def get_kinematics_from_int(int_variables, constants:dict, return_p_in=False, return_spherical=False, epsilon:float=1e-3, err_map=err_map, do_check=True, dEmax:float=1., dpmax:float=1., nwa=True):
    """Reference implementation of kinematics solver
    Available in C++ and callable from Python in lib.py
    as calc_kinematics_from_int

    Args:
        int_variables (_type_): _description_
        constants (dict): _description_
        return_p_in (bool, optional): _description_. Defaults to False.
        return_spherical (bool, optional): _description_. Defaults to False.
        epsilon (float, optional): _description_. Defaults to 1e-3.
        err_map (_type_, optional): _description_. Defaults to err_map.
        do_check (bool, optional): _description_. Defaults to True.
        dEmax (float, optional): _description_. Defaults to 1..
        dpmax (float, optional): _description_. Defaults to 1..
        nwa (bool, optional): _description_. Defaults to True.

    Returns:
        _type_: _description_
    """
    
    mH2 = Thb1 = Phb1 = Rhb1 = Thb1b = Phb1b = Rhb2 = Thb2 = 0
    if nwa:
        mH2 = 125.**2
        Thb1, Phb1, Rhb1, Thb1b, Phb1b, Rhb2, Thb2 = int_variables
    else:
        mH2, Thb1, Phb1, Rhb1, Thb1b, Phb1b, Rhb2, Thb2 = int_variables
    
    #print("Thb2:", Thb2)
    #print("Rhb2:", Rhb2)
    
    b1E = sqrt(constants["m_b"]**2 + Rhb1**2)
    
    b1e = unit_vec(Thb1, Phb1)
    b1p = Rhb1*b1e
    
    # Calculate Rhb1b from mH2, Rhb1 and angles
    b1be = unit_vec(Thb1b, Phb1b)
    dot = np.dot(b1e, b1be)
    
    dmH2 = (mH2-2*constants["m_b"]**2)/2
    d = (dmH2*Rhb1*dot)/(Rhb1**2*dot**2 - b1E**2)
    
    arg_Rhb1b_sqrt = ((b1E**2)*(constants["m_b"]**2) - dmH2**2)/((Rhb1**2)*(dot**2) - b1E**2) + d**2
    if arg_Rhb1b_sqrt < 0:
        return err_map[1]
    
    #Rhb1b2 = -d -sqrt(arg_Rhb1b_sqrt) # negative solution
    Rhb1b = -d +sqrt(arg_Rhb1b_sqrt)
    if Rhb1b < 0:
        return err_map[2]
    
    b1bp1 = Rhb1b*b1be
    #b1bp2 = Rhb1b2*b1be
    
    b1bE1 = sqrt(constants["m_b"]**2 + Rhb1b**2)
    #print("b1E", b1E)
    #print("b1bE", b1bE1)
    
    # Calculate pH and from that pB2
    pB2 = np.array([
        constants["system_E"] -(b1E+b1bE1),
        constants["system_px"] -b1p[0] -b1bp1[0],
        constants["system_py"] -b1p[1] -b1bp1[1],
        constants["system_pz"] -b1p[2] -b1bp1[2],
    ])
    #print("pB2", pB2)
    
    # Calculate Rhb2b
    b2E = sqrt(constants["m_b"]**2 + Rhb2**2)
    #print("b2E", b2E)
    #print("b2bE: pB2[0]-b2E", pB2[0], "-", b2E, "=", pB2[0]-b2E)
    
    b2bE = pB2[0] - b2E
    if b2bE < 0:
        return err_map[5]
    #print("b2bE", b2bE)
    
    arg_Rhb2b_sqrt = b2bE**2 - constants["m_b"]**2
    if arg_Rhb2b_sqrt < 0:
        return err_map[6]
    Rhb2b = sqrt(arg_Rhb2b_sqrt)
    #print("Rhb2b", Rhb2b)
    
    # Calculate remaining variables, i.e. Thb2b, Phb2b and Phb2
    arg_Thb2b_acos = 1/Rhb2b*(pB2[3] - Rhb2*cos(Thb2))
    if abs(arg_Thb2b_acos) > 1:
        return err_map[7]
    Thb2b = acos(arg_Thb2b_acos)
    
    a = pB2[1]
    b = pB2[2]
    c = (-(Rhb2*sin(Thb2))**2 + (Rhb2b*sin(Thb2b))**2 + a**2 + b**2)/(2*Rhb2b*sin(Thb2b))
    
    arg_Phb2b_atan_sqrt = a**2 + b**2 - c**2
    if arg_Phb2b_atan_sqrt < 0:
        if -arg_Phb2b_atan_sqrt < epsilon:
            arg_Phb2b_atan_sqrt = 0
        else:
            return err_map[8]
    
    # Two solutions to calculate Phb2b with atan; try out both, use the one that is closest to target   
    Phb2b_1 = 2*atan((b+ sqrt(arg_Phb2b_atan_sqrt))/(a+c))
    Phb2b_2 = 2*atan((b- sqrt(arg_Phb2b_atan_sqrt))/(a+c))
    
    arg_Phb2_1_acos = (a - Rhb2b*sin(Thb2b)*cos(Phb2b_1))/(Rhb2*sin(Thb2))
    arg_Phb2_2_acos = (a - Rhb2b*sin(Thb2b)*cos(Phb2b_2))/(Rhb2*sin(Thb2))
    
    if abs(arg_Phb2_1_acos) and abs(arg_Phb2_2_acos) > 1:
        return err_map[9]
    
    # Check which solution is better    
    Phb2 = 0
    Phb2b = 0
    if abs(arg_Phb2_1_acos) > 1:
        Phb2b = Phb2b_2
        Phb2 = acos(arg_Phb2_2_acos)
        
        if (pB2[2] - Rhb2b*sin(Thb2b)*sin(Phb2b))/(sin(Thb2)) < 0:
            Phb2 = -Phb2
    elif abs(arg_Phb2_2_acos) > 1:
        Phb2b = Phb2b_1
        Phb2 = acos(arg_Phb2_1_acos)
        
        if (pB2[2] - Rhb2b*sin(Thb2b)*sin(Phb2b))/(sin(Thb2)) < 0:
            Phb2 = -Phb2
    else:
        Phb2_1 = acos(arg_Phb2_1_acos)
        Phb2_2 = acos(arg_Phb2_2_acos)
        
        #print("Sign Hint")
        #print((pB2[1] - Rhb2b*sin(Thb2b)*cos(Phb2b_1))/(sin(Thb2)))
        #print((pB2[1] - Rhb2b*sin(Thb2b)*cos(Phb2b_2))/(sin(Thb2)))
        
        if (pB2[2] - Rhb2b*sin(Thb2b)*sin(Phb2b_1))/(sin(Thb2)) < 0:
            Phb2_1 = -Phb2_1
            
        if (pB2[2] - Rhb2b*sin(Thb2b)*sin(Phb2b_2))/(sin(Thb2)) < 0:
            Phb2_2 = -Phb2_2
        
        b2_1p = Rhb2*unit_vec(Thb2, Phb2_1)
        b2_2p = Rhb2*unit_vec(Thb2, Phb2_2)
        
        b2b_1p = Rhb2b*unit_vec(Thb2b, Phb2b_1)
        b2b_2p = Rhb2b*unit_vec(Thb2b, Phb2b_2)
        
        B1comp_1 = b2_1p + b2b_1p
        B1comp_2 = b2_2p + b2b_2p
        
        B1check_1 = pB2[1:4] - B1comp_1
        B1check_2 = pB2[1:4] - B1comp_2
        
        if np.dot(B1check_1, B1check_1) < np.dot(B1check_2, B1check_2):
            Phb2 = Phb2_1
            Phb2b = Phb2b_1
        else:
            Phb2 = Phb2_2
            Phb2b = Phb2b_2

    # Finish
    b2e = unit_vec(Thb2, Phb2)
    b2p = Rhb2*b2e
    
    b2be = unit_vec(Thb2b, Phb2b)
    b2bp = Rhb2b*b2be
    
    result = np.array([
        (b1E, *b1p),
        (b1bE1, *b1bp1),
        (b2E, *b2p),
        (b2bE, *b2bp)
    ])
    
    if do_check:
        check = np.array([
            constants["system_E"] - np.sum(result.T[0]),
            constants["system_px"] - np.sum(result.T[1]),
            constants["system_py"] - np.sum(result.T[2]),
            constants["system_pz"] - np.sum(result.T[3])
        ])
        
        dE = check[0]
        dp = sqrt(check[1]**2 + check[2]**2 + check[3]**2)
        
        if dE > dEmax or dp > dpmax:
            return err_map[10]
    
    if return_spherical:
        return (
            np.array([
                Rhb1, Thb1, Phb1,
                Rhb1b, Thb1b, Phb1b,
                Rhb2, Thb2, Phb2,
                Rhb2b, Thb2b, Phb2b
            ]),
            result
        )
    
    if not return_p_in:
        return result
    else:
        return (
            np.array([
                constants["system_E"],
                constants["system_px"],
                constants["system_py"],
                constants["system_pz"]
            ]),
            result
        )

def jac(params, constants):
    # Jacobian for transformation from
    # {Rhb1b, Rhb2b, Thb2b, Phb2, Phb2b, Rhb1, Rhb2, Thb1, Thb1b, Phb1, Phb1b, Thb2} to
    #              {PE, Px, Py, Pz, mH2, Rhb1, Rhb2, Thb1, Thb1b, Phb1, Phb1b, Thb2}
    
    (Rhb1, Thb1, Phb1,
    Rhb1b, Thb1b, Phb1b,
    Rhb2, Thb2, Phb2,
    Rhb2b, Thb2b, Phb2b) = params
    
    # Theta    
    mb = constants["m_b"]
    
    return abs(
        (1/sqrt(mb**2 + Rhb2b**2)) * (
            Rhb2*(Rhb2b**3)*sin(Phb2 - Phb2b)* 2*(
                + sqrt(mb**2 + Rhb1**2)*Rhb1b/sqrt(mb**2 + Rhb1b**2)
                - Rhb1*(
                      cos(Thb1)*cos(Thb1b)
                    + cos(Phb1-Phb1b)*sin(Thb1)*sin(Thb1b)
                )
            )*sin(Thb2)*sin(Thb2b)**2
        )
    )

def construct_integrand_bf(data, event_idx, constants):
    from analysis.cffi.mg5.lib import calc_zhh
    
    muon_kin = prepare_int(data, event_idx, constants)
    ((mu1E, mu1p), (mu2E, mu2p)) = muon_kin
    
    evt_constants = get_evt_constants(data, event_idx, constants)
    
    reco_energies, reco_momenta, reco_angles = get_kinematics_tf(data, event_idx)
    
    @vegas.batchintegrand
    def integrand(vars):
        results = np.ones(len(vars), dtype=np.float64)
        valid_idx = []
        valid_kinematics = []
        found = 0
        
        for i in range(len(vars)):
            x = vars[i]
        
            result = get_kinematics_from_int(x, evt_constants, return_spherical=True)

            if isinstance(result, tuple):
                # print("Found valid phase space point")
                
                spherical, kin = result
                
                # Valid kinematics -> evaluate matrix element and jacobian
                full_kin = np.array([
                    [mu1E, *mu1p],
                    [mu2E, *mu2p],
                    kin[0],
                    kin[1],
                    kin[2],
                    kin[3]
                ])
                
                valid_kinematics.append(full_kin)
                valid_idx.append(i)
                
                jacobian = jac(spherical, constants)
                
                # Transfer functions and parton level quantities
                spherical = np.reshape(spherical, (4,3))
                rhos = spherical.T[0]
                
                energies = np.sqrt([
                    constants["m_b"]**2 + rhos**2
                ])[0]
                
                thetas = spherical.T[1]
                phis = spherical.T[2]
                
                transfer = 1
                for j in range(4):
                    transfer = transfer * tf_E(reco_energies[j+2] - energies[j]) * tf_Th(reco_angles[j+2][0] - thetas[j]) * tf_Ph(reco_angles[j+2][1] - phis[j])
                
                results[i] = jacobian*transfer
                found = found+1
                
            else:
                # No valid kinematics -> reject
                results[i] = 0
        
        print(f"Found PS points {found}/{len(vars)} ({(100*found/len(vars)):.2f} %)")
        
        # Calculate matrix elements in batch mode        
        valid_kinematics = np.array(valid_kinematics).flatten().tolist()
        melements = calc_zhh(valid_kinematics)
        
        np.put(results, valid_idx, results[valid_idx]*melements)
        
        if False:
            print(len(valid_kinematics), len(valid_idx))
        
            exp_from = -33
            exp_to = -18
            
            plt.hist(results[results != 0], bins=np.logspace(exp_from, exp_to, 64))
            #plt.yscale('log')
            plt.xscale('log')
            plt.xlim((10**exp_from, 10**exp_to))
            #plt.ylim((0., 0.1))true_lep1_pxtrue_lep1_pxtrue_lep1_pxtrue_lep1_pxtrue_lep1_pxtrue_lep1_px
            plt.show()
                
        return results
    
    return integrand


def int_bf(data, event_idx:int, evt_constants:dict, size:int=100000):
    """MEM integration using the Python integrand

    Args:
        data (_type_): _description_
        event_idx (_type_): _description_
        evt_constants (dict): _description_
        size (int, optional): _description_. Defaults to 100000.
    """
    
    muon_kin = prepare_int(data, event_idx, evt_constants, True)
    ((mu1E, mu1p), (mu2E, mu2p)) = muon_kin
    
    # phi: from -pi to pi
    # theta: from 0 to pi
    
    bounds = [
        [124.9**2, 125.1**2],
        [0, pi],
        
        [-pi, pi],
        [0, 200],
        
        [0, pi],
        [-pi, pi],
        
        [0, 200],
        [0, pi]
    ]
    
    map = vegas.AdaptiveMap(bounds)
    x = np.stack([
        np.ones(size)*125**2, # mB1pow2
        Thb_prior(size), # Thb1
        
        Phb_prior(size), # Phb1
        Rhb_prior(size), # Rhb1
        
        Thb_prior(size), # Thb1b
        Phb_prior(size), # Phb1b
        
        Rhb_prior(size), # Rhb2
        Thb_prior(size), # Thb2
    ]).T
    
    f = construct_integrand_bf(data, 0, evt_constants)
    
    map.adapt_to_samples(x, f(x), nitn=5)
    
    integ = vegas.Integrator(map, alpha=0.1)
    result = integ(f, nitn=10, neval=1000000, save=f'zhh_{event_idx}.pkl')
    
    print(result.summary())
    print('result = %s    Q = %.2f' % (result, result.Q))
    
    return result


def int_test_bf(data, event_idx:int, constants:dict, samples_per_dim:int=4, use_reco=True):   
    muon_kin = prepare_int(data, event_idx, constants, use_reco)
    ((mu1E, mu1p), (mu2E, mu2p)) = muon_kin
    
    boundaries = [
        (124.9**2, 125.1**2), # mH2
        (0., pi), # Thb1
        
        (-pi, pi), # Phb1
        (0., 200.), # Rhb1
        
        (0., pi), # Thb1b
        (-pi, pi), # Phb1b
        
        (0., 200.), # Rhb2
        (0., pi) # Thb2
    ]

    samples = []
    for boundary in boundaries:
        samples.append(np.linspace(boundary[0], boundary[1], samples_per_dim))
        
    print(samples)
    
    evt_constants = constants | {
        "sqrt_s": constants["sqrt_s"] -mu1E -mu2E,
        "system_px": constants["system_px"] -mu1p[0] -mu2p[0],
        "system_py": constants["system_py"] -mu1p[1] -mu2p[1],
        "system_pz": constants["system_pz"] -mu1p[2] -mu2p[2],
    }
    
    values = np.zeros(8)
    
    valid = 0
    fail = 0
    
    kins = []
    errors = []
    inputs = []
    for i in product(range(samples_per_dim), repeat=8):
        for j in range(8):
            values[j] = samples[j][i[j]]
        
        # evaluate
        # arg order: (mu1E, mu1p), (mu2E, mu2p), mB1pow2, Thb1, Phb1, Rhb1, Thb1b, Phb1b, Rhb2, Thb2 = int_variables
        try:
            result = get_kinematics_from_int(
                values, evt_constants, return_p_in=True
            )
            
            if isinstance(result, tuple):
                pIn, kin = result
                
                valid += 1
                kins.append(kin)
                inputs.append(values.copy())
            else:
                fail += 1
                errors.append(result)
            
        except Exception:
            #print(f"Exception occured with params=", values, "at i=", i)
            fail += 1
            errors.append(-99)
    
    print("Done")
    
    return valid, fail, kins, errors, inputs