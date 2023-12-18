import click
from analysis import logger

@click.group(name="mem")
def cli_mem():
    """Commands for MEM integration"""
    pass

@cli_mem.command()
@click.argument("event")
@click.argument("dst")
@click.option("--me_type", default="1", help="0 for MG5; 1 for Physsim", type=int)
@click.option("--src", default="/nfs/dust/ilc/user/bliewert/fullflow_v3/comparison/cache/comparison_reco_zhh_zzh.npy", help="Numpy file with reco kinematics")
def integrate(event:int, dst:str, me_type:int, src:str):
    """Does the MEM integration for the given event for signal and background hypothesis and saves the results in dst

    Args:
        event (int): _description_
        dst (str): _description_
        me_type (int): 
        src (str): _description_
    """
    from analysis.mem import int_bf_v2
    from analysis.import_data import import_true_reco
    import numpy as np
    import pandas as pd
    
    event = int(event)
    
    reco = import_true_reco(src_file=src, normalize=False)
    
    event_idx = np.where(reco["event"] == event)[0]
    if len(event_idx) == 0:
        raise Exception(f"Event [{event}] could not be found")
    
    event_idx = event_idx[0]
    
    f = open(dst, "a")
    f.write(f"Event [{event}] IDX [{event_idx}] \n")
    f.close()
    
    logger.info("Starting ZHH")
    
    neval = 10000000
    nitn = 10
    
    #if me_type == 0:
    #    neval = 10*neval
    #    nitn = 4
    
    res_sig = int_bf_v2(reco, event_idx, mode=1, neval=neval, precond_size=2000000, nitn=nitn, me_type=me_type)
    
    f = open(dst, "a")
    f.write("ZHH: " + str(res_sig) + "\n")
    f.close()
    
    logger.info("Finished ZHH")
    
    
    logger.info("Starting ZZH")
    
    neval = 500000
    nitn = 10
    
    res_bkg = int_bf_v2(reco, event_idx, mode=0, neval=neval, precond_size=2000000, nitn=nitn, me_type=me_type)
    
    f = open(dst, "a")
    f.write("ZZH: " + str(res_bkg) + "\n")
    f.close()
    
    logger.info("Finished ZZH")


if __name__ == '__main__':
    cli_mem()