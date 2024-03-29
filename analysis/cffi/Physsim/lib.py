from analysis.cffi.Physsim.CalculateME import lib as Physsim
from typing import List

def inv_vec(vec:List[float]):
    return [vec[1:4], vec[0]]

def vec_sum(vec1:List[float], vec2:List[float]):
    return [
        vec1[0] + vec2[0],
        vec1[1] + vec2[1],
        vec1[2] + vec2[2],
        vec1[3] + vec2[3],
    ]

def calculate_me(reco_kin:List[float],
                 mode:int,
                 )->float:
    """Calculate the ZHH/ZZH matrix element with Physsim
    Physsim calculates sqrt_s from the kinematics

    Args:
        reco_kin (List[float]): _description_
        mode (int): 1 ZHH, 0 ZZH

    Raises:
        Exception: _description_

    Returns:
        _type_: _description_
    """
    
    if not isinstance(reco_kin, list):
        raise Exception('reco_kin must be list')
    
    if mode != 1 and mode != 0:
        raise Exception('Invalid mode')
    
    if mode == 0:
        kin = [
            reco_kin[0 ]                , reco_kin[1 ]              , reco_kin[2 ]              , reco_kin[3 ],
            reco_kin[4 ]                , reco_kin[5 ]              , reco_kin[6 ]              , reco_kin[7 ],
            reco_kin[8 ]                , reco_kin[9 ]              , reco_kin[10]              , reco_kin[11],
            reco_kin[12]                , reco_kin[13]              , reco_kin[14]              , reco_kin[15],
            reco_kin[16]+reco_kin[20]   , reco_kin[17]+reco_kin[21] , reco_kin[18]+reco_kin[22] , reco_kin[19]+reco_kin[23],
        ]
        res = Physsim.calc_zzh_single(kin)
    else:
        kin = [
            reco_kin[0 ]                , reco_kin[1 ]              , reco_kin[2 ]              , reco_kin[3 ],
            reco_kin[4 ]                , reco_kin[5 ]              , reco_kin[6 ]              , reco_kin[7 ],
            reco_kin[8 ]+reco_kin[12]   , reco_kin[9 ]+reco_kin[13] , reco_kin[10]+reco_kin[14] , reco_kin[11]+reco_kin[15],
            reco_kin[16]+reco_kin[20]   , reco_kin[17]+reco_kin[21] , reco_kin[18]+reco_kin[22] , reco_kin[19]+reco_kin[23],
        ]
        res = Physsim.calc_zhh_single(kin)
    
    return res