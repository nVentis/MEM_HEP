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
@click.option("--src", default="/nfs/dust/ilc/user/bliewert/fullflow_v3/comparison/npy/reco/compare_reco_with_mg5.npy", help="Numpy file with reco kinematics")
def integrate(event:int, dst:str, me_type:int, src:str):
    """Does the MEM integration for the given event for signal and background hypothesis and saves the results in dst

    Args:
        event (int): _description_
        dst (str): _description_
        me_type (int): 
        src (str): _description_
    """
    from analysis.mem import int_bf_v2
    import numpy as np
    import pandas as pd
    
    event = int(event)
    
    reco = pd.DataFrame(np.load(src, allow_pickle=True))
    
    f = open(dst, "a")
    f.write(f"Event [{event}] \n")
    f.close()
    
    logger.info("Starting ZHH")
    
    neval = 6000000
    nitn = 10
    
    #if me_type == 0:
    #    neval = 10*neval
    #    nitn = 4
    
    res_sig = int_bf_v2(reco, event, mode=1, neval=neval, precond_size=2000000, nitn=nitn, me_type=me_type)
    
    f = open(dst, "a")
    f.write("ZHH: " + str(res_sig) + "\n")
    f.close()
    
    logger.info("Finished ZHH")
    
    
    logger.info("Starting ZZH")
    
    neval = 1000000
    nitn = 10
    
    res_bkg = int_bf_v2(reco, event, mode=0, neval=neval, precond_size=2000000, nitn=nitn, me_type=me_type)
    
    f = open(dst, "a")
    f.write("ZZH: " + str(res_bkg) + "\n")
    f.close()
    
    logger.info("Finished ZZH")


if __name__ == '__main__':
    cli_mem()