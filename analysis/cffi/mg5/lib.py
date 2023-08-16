from analysis.cffi.mg5.CalcMEZHH import lib as zhh
from analysis.cffi.mg5.CalcMEZZH import lib as zzh
from typing import List

def calc_zhh(momenta: List[float],
                    helicities = [0,-1,1,1],
                    energy:float = 500.,
                    param_card:str = "/afs/desy.de/user/b/bliewert/public/MarlinWorkdirs/MEM_HEP/analysis/cffi/mg5/mg5/Cards/param_card.dat") -> float:
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
                    helicities = [0,-1,1,1],
                    energy:float = 500.,
                    param_card:str = "/afs/desy.de/user/b/bliewert/public/MarlinWorkdirs/MEM_HEP/analysis/cffi/mg5/mg5/Cards/param_card.dat") -> float:
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