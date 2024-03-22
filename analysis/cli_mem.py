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
@click.option("--sampling", type=click.Choice(['vegas', 'nis'], case_sensitive=False))
@click.option("--src", default="/nfs/dust/ilc/user/bliewert/fullflow_v3/comparison/cache/comparison_reco_zhh_zzh.npy", help="Numpy file with reco kinematics")
def integrate(event:int, dst:str, me_type:int, sampling_strategy:str, src:str):
    """Does the MEM integration for the given event for signal and background hypothesis and saves the results in dst

    Args:
        event (int): _description_
        dst (str): _description_
        me_type (int): 
        src (str): _description_
    """
    from analysis.mem import mem_integrate
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
    
    res_sig = mem_integrate(reco, event_idx, mode=1, neval=neval, precond_size=2000000, nitn=nitn, me_type=me_type)
    logger.info(f"Finished ZHH: [{str(res_sig)}]")
    
    
    logger.info("Starting ZZH")
    
    neval = 500000
    nitn = 10
    
    res_bkg = mem_integrate(reco, event_idx, mode=0, neval=neval, precond_size=2000000, nitn=nitn, me_type=me_type)
    
    report = f"""
ZHH: {str(res_sig)}
ZZH: {str(res_bkg)}             """

    f = open(dst, "a")
    f.write(report)
    f.close()
    
    logger.info(f"Finished ZZH: [{str(res_bkg)}]")


if __name__ == '__main__':
    cli_mem()