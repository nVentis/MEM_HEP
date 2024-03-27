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
@click.option("--save_npy", default="1", type=int)
@click.option("--src", default="/nfs/dust/ilc/user/bliewert/fullflow_v3/comparison/cache/comparison_reco_zhh_zzh.npy", help="Numpy file with reco kinematics")
def integrate(event:int, dst:str, me_type:int, sampling:str, save_npy:int, src:str):
    """Does the MEM integration for the given event for signal and background hypothesis and saves the results in dst
    Careful naming: event (input) is the index in the data frame, event_idx is the generator event id

    Args:
        event (int): _description_
        dst (str): _description_
        me_type (int): 
        src (str): _description_
    """
    from analysis.mem import mem_integrate
    from analysis.import_data import import_true_reco
    from analysis.calc import get_kinematics
    import itertools
    import numpy as np
    import pandas as pd
    from os.path import dirname, join
    
    event = int(event)
    me_type = int(me_type)
    save_npy = bool(save_npy)
    int_type = 0 if sampling == "vegas" else 1
    
    reco = import_true_reco(src_file=src, normalize=False)
    event_idx = int(reco.iloc[event]['event'])
    
    perms_zhh = [
        [1,2,3,4],
        [1,3,2,4],
        [1,4,2,3],
    ]
    perms_zzh = list(itertools.permutations([1, 2, 3, 4]))
    
    ########################################################
    # Helper function and logging
    def extract_values(result, int_type:int):
        """Extract the values to save for the different integration methods

        Args:
            result (_type_): _description_
            int_type (int): _description_

        Raises:
            Exception: _description_

        Returns:
            _type_: _description_
        """
        if int_type == 0:
            mean, stddev = result.mean, result.sdev
            
            return ((mean, stddev), result.itn_results)
        elif int_type == 1:
            mean, stddev = result['result'], result['tot_uncert']
            
            return ((mean, stddev), result)
        else:
            raise Exception('Invalid state')
    
    f = open(dst, "a")
    f.write(f"Event [{event}] IDX [{event_idx}] \n")
    f.close()
    
    precond_size = 2000000
    
    ########################################################
    # ZHH
    logger.info("Starting ZHH")
    
    neval = 10000000
    nitn = 10
    
    #if me_type == 0:
    #    neval = 10*neval
    #    nitn = 4
    
    results_sig = { 'means': [], 'sigmas': [] }
    results_bkg = { 'means': [], 'sigmas': [] }
    
    for i in range(len(perms_zhh)):
        perm = perms_zhh[i]
        logger.info(f'Perm [{i}]: [{" ".join(str(x) for x in perm)}]')
        
        reco_kin = get_kinematics(reco, False, i=event, perm=perm)
        
        result = mem_integrate(reco_kin=reco_kin, mode=1, neval=neval, precond_size=precond_size, nitn=nitn, me_type=me_type, int_type=int_type, file_suffix=str(i))
        (mean, sdev) , save_result = extract_values(result, int_type=int_type)
        
        results_sig['means'].append(mean)
        results_sig['sigmas'].append(sdev)
        logger.info(f'Result: {str(result)}')
        
        if save_npy:
            np.save(join(dirname('./result.txt'), f'sig_p{str(i)}.npy'), np.array(save_result), allow_pickle=True)
    
    tot_sig = np.mean(results_sig['means'])
    
    logger.info(f"Finished ZHH: [{str(tot_sig)}]")
    
    ########################################################
    # ZZH
    logger.info("Starting ZZH")
    
    neval = 500000
    nitn = 10
    
    for i in range(len(perms_zzh)):
        perm = perms_zzh[i]
        logger.info(f'Perm [{i}]: [{" ".join(str(x) for x in perm)}]')
        
        reco_kin = get_kinematics(reco, False, i=event, perm=perm)
        
        result = mem_integrate(reco_kin=reco_kin, mode=0, neval=neval, precond_size=precond_size, nitn=nitn, me_type=me_type, int_type=int_type, file_suffix=str(i))
        (mean, sdev) , save_result = extract_values(result, int_type=int_type)
        
        results_bkg['means'].append(mean)
        results_bkg['sigmas'].append(sdev)
        logger.info(f'Result: {str(result)}')
        
        if save_npy:
            np.save(join(dirname('./result.txt'), f'bkg_p{str(i)}.npy'), np.array(save_result), allow_pickle=True)
    
    tot_bkg = np.mean(results_bkg['means'])
    logger.info(f"Finished ZZH: [{str(tot_bkg)}]")
    
    ########################################################
    # Save results
    
    if save_npy:
        np.save(join(dirname('./result.txt'), 'summary.npy'), {
            'sig': results_sig,
            'bkg': results_bkg
        }, allow_pickle=True)
    
    report = f"""
ZHH: [{str(tot_sig)}]
ZZH: [{str(tot_bkg)}]
"""

    f = open(dst, "a")
    f.write(report)
    f.close()
    
    logger.info(f'Done')


if __name__ == '__main__':
    cli_mem()