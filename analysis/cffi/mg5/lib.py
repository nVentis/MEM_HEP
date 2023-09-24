from analysis.cffi.mg5.CalcMEZHH import lib as zhh
from analysis.cffi.mg5.CalcMEZZH import lib as zzh
from analysis.cffi.mg5.compiled_with import lib_options
from typing import List

def mc_batch(reco_kin:List[float],
                 int_variables:List[float],
                 mode:int,
                 helicities:List[int] = [0,-1,1,1],
                 energy:float = 500.,
                 param_card:str = "/afs/desy.de/user/b/bliewert/public/MarlinWorkdirs/MEM_HEP/analysis/cffi/mg5/mg5/Cards/param_card.dat",
                 tf_E_params=None,
                 tf_Th_params=None,
                 tf_Ph_params=None,
                 NWA:bool=lib_options["NWA"]):
    """Monte Carlo integrand

    Args:
        reco_kin (_type_): expects kinematics in (E,px,py,pz) form for (mu-,mu+,b,b,b,b)
        int_variables (_type_): integration variables
        mode (int): 1 for zhh 0 for zzh
        helicities (list, optional): _description_. Defaults to [0,-1,1,1].
        energy (float, optional): _description_. Defaults to 500..
        param_card (str, optional): _description_. Defaults to "/afs/desy.de/user/b/bliewert/public/MarlinWorkdirs/MEM_HEP/analysis/cffi/mg5/mg5/Cards/param_card.dat".
        tf_E_params,tf_Th_params,tf_Ph_params: (x0,gamma) parameters for Lorentzian transfer functions. if None, defaults are used. Defaults to None.
        NWA (bool, optional):  
    """
    
    n_elements = int(len(int_variables)/(7 if NWA else 8))
    
    lib = (zhh if mode else zzh)
    result = lib.calc_mc_batch(str.encode(param_card), energy, helicities, int(len(helicities)/2), reco_kin, int_variables, n_elements)
    
    res_list = [result[i] for i in range(n_elements)]
    
    lib.free(result)
    
    return res_list
    
def calc_kinematics_from_int(
                int_variables:List[float],
                evt_constants:List[float],
                helicities:List[int] = [0,-1,1,1],
                energy:float = 500.,
                param_card:str = "/afs/desy.de/user/b/bliewert/public/MarlinWorkdirs/MEM_HEP/analysis/cffi/mg5/mg5/Cards/param_card.dat",
                NWA:bool=lib_options["NWA"]):
    """Monte Carlo integrand

    Args:
        int_variables (List[float]): integration variables (mH2, Thb1, Phb1, Rhb1, Thb1b, Phb1b, Rhb2, Thb2); without mH2 if NWA is True
        evt_constants (List[float]): mb, epsilon, dEmax, dpmax, sqrt_s, system_E, system_px, system_py, system_pz
        helicities (list, optional): _description_. Defaults to [0,-1,1,1].
        energy (float, optional): _description_. Defaults to 500..
        param_card (str, optional): _description_. Defaults to "/afs/desy.de/user/b/bliewert/public/MarlinWorkdirs/MEM_HEP/analysis/cffi/mg5/mg5/Cards/param_card.dat".
        tf_E_params,tf_Th_params,tf_Ph_params: (x0,gamma) parameters for Lorentzian transfer functions. if None, defaults are used. Defaults to None.
        NWA (bool): Use Narrow Width approximation; mH2 implicitly set to 125.**2 in C++ code
    """
    
    if len(int_variables) != (7 if NWA else 8):
        raise Exception("Invalid format of int_variables")
    
    if len(evt_constants) != 9:
        raise Exception("Invalid format of evt_constants")
    
    return zhh.calc_kinematics_from_int(str.encode(param_card), evt_constants, helicities, int(len(helicities)/2), *int_variables)
    
def calc_zhh(momenta: List[float],
                    helicities:List[int] = [0,-1,1,1],
                    energy:float = 500.,
                    param_card:str = "/afs/desy.de/user/b/bliewert/public/MarlinWorkdirs/MEM_HEP/analysis/cffi/mg5/mg5/Cards/param_card.dat") -> List[float]:
    """Calculates the e-e+ -> ZHH -> mu- mu+ b bbar b bbar squared matrix element using mg5 for m events given momenta and helicity states

    Args:
        momenta (List[float]): flat list length m of (E,Px,Py,Pz) for l particles (lep1, lep2, b1, bbar1, b2, bbar2), concatenated to length m=6*l
        helicities (list, optional): flat list of 2n entries of where (particle position)@2n and (helicity +/-1)@(2n+1). Defaults to [0, -1, 1, 1].
        energy (float, optional): COM energy; kinematics assumed to be aligned along Z 
        param_card (str, optional): physics parameter file (widths, masses etc). Defaults to "/afs/desy.de/user/b/bliewert/public/MarlinWorkdirs/MEM_HEP/analysis/cffi/mg5/mg5/Cards/param_card.dat".

    Returns:
        float: the squared matrix element
    """
    
    n_elements = int(len(momenta)/24)
    
    result = zhh.calc(str.encode(param_card), energy, helicities, int(len(helicities)/2), momenta, n_elements)
    res_list = [result[i] for i in range(n_elements)]
    
    zhh.free(result)
    
    return res_list

def calc_zzh(momenta: List[float],
                    helicities:List[int] = [0,-1,1,1],
                    energy:float = 500.,
                    param_card:str = "/afs/desy.de/user/b/bliewert/public/MarlinWorkdirs/MEM_HEP/analysis/cffi/mg5/mg5/Cards/param_card.dat") -> List[float]:
    """Calculates the e-e+ -> ZZH -> mu- mu+ b bbar b bbar squared matrix element using mg5 for m events given momenta and helicity states

    Args:
        momenta (List[float]): flat list length m of (E,Px,Py,Pz) for l particles (lep1, lep2, b1, bbar1, b2, bbar2), concatenated to length m=6*l
        helicities (list, optional): flat list of 2n entries of where (particle position)@2n and (helicity +/-1)@(2n+1). Defaults to [0, -1, 1, 1].
        energy (float, optional): COM energy; kinematics assumed to be aligned along Z 
        param_card (str, optional): physics parameter file (widths, masses etc). Defaults to "/afs/desy.de/user/b/bliewert/public/MarlinWorkdirs/MEM_HEP/analysis/cffi/mg5/mg5/Cards/param_card.dat".

    Returns:
        float: the squared matrix element
    """
    
    n_elements = int(len(momenta)/24)
    
    result = zzh.calc(str.encode(param_card), energy, helicities, int(len(helicities)/2), momenta, n_elements)
    res_list = [result[i] for i in range(n_elements)]
    
    zzh.free(result)
    
    return res_list

if __name__ == "__main__":
    print("Test")
    
    fm = [
        20., 1, 2, 3,
        25, 2, 3, 4,
        100, 3, 4, 5,
        120, 4, 5, 6,
        150, 5, 6, 7,
        160, 6, 7, 8, 
    ]

    fmB = [
        25., 1, 2, 3,
        30, 2, 3, 4,
        105, 3, 4, 5,
        125, 4, 5, 6,
        155, 5, 6, 7,
        165, 6, 7, 8, 
    ]

    fm2 = fm + fmB

    from time import time

    t_s = time()
    for i in range(1000):
        calc_zhh(fm)

    t_s1 = time() - t_s

    fmt = fm
    for i in range(1000):
        fmt = fmt + fm

    t_s = time()

    calc_zhh(fmt)

    t_s2 = time() - t_s

    print(t_s1, t_s2)

    #zhh.calc(b"/afs/desy.de/user/b/bliewert/public/MarlinWorkdirs/MEM_HEP/analysis/cffi/mg5/mg5/Cards/param_card.dat", 2, [0,-1,1,1])
    #zzh.calc(b"/afs/desy.de/user/b/bliewert/public/MarlinWorkdirs/MEM_HEP/analysis/cffi/mg5/mg5/Cards/param_card.dat", 2, [0,-1,1,1])