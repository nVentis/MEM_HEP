import click

import sys
import os
import shutil
import logging
logger = logging.getLogger('mem_hep')

# Expects two arguments:
# 1. comparison_root file
# 2. comparison_out directory: Where to store the converted data and plots
# 3. optional: production plot mode; if true, uses ROOT, otherwise uses matplotlib

from analysis.convert import convert_file
from analysis.import_data import import_data, filter_data, combine_columns, split_true_zhh_zzh
from analysis.calc import calc_nll_llr_dtf_delta, calc_nll_llr_dtf_dbgauss
from analysis.plot_routines import plot_summary, plot_nll, plot_llr

def me(src_file:str, dst:str, name:str, convert:bool, plot:bool, calc:bool = True, dtf:str = "delta"):
    """Attempts separating ZHH vs ZZH using matrix elements"""
    logger.info("Analyzing {}".format(name))
    
    root_file = os.path.basename(src_file)
    cnv_file = dst if os.path.splitext(root_file)[1] == "npy" else os.path.join(dst, os.path.splitext(root_file)[0] + ".npy")
    dst_dir = os.path.dirname(cnv_file)
    
    if convert:
        from analysis.convert import convert_file
        
        convert_file(src_file, "dataTree", cnv_file)
        logger.info("Converted file {} to {}".format(src_file, cnv_file))
    else:
        logger.info("Skipped conversion")
        
    if plot or calc:
        raw = import_data(cnv_file)

        if plot:
            plots_dir = os.path.join(dst_dir, "plots")
            for dir in [plots_dir]:
                if os.path.isdir(dir):
                    shutil.rmtree(dir, ignore_errors=True)
                    logger.info("Clearing directory-tree {}".format(dir))
                    
                os.makedirs(dir)

            # Show occurences of error_code        
            summary_path = os.path.join(dst_dir, "plots", "summary.png")
            plot_summary(raw, name, summary_path)
            logger.info("Saved error summary at {}".format(summary_path))
            
            # Now plot nll and llr
            # Filter to only contain entries without errors (error_code = 0) and with matrix elements > 0 for ZHH and ZZH (zhh_sigmalr and zzh_sigmalr)
            filtered = filter_data(raw)
            
            # A Delta-distributions as transfer functions
            
            # First, calculate the log likelihood ratio in case of delta distributions as transfer functions,
            # i.e. assuming measured properties are parton-level-properties; useful only for MCParticle/TrueJet assumption
            # Breaks down when real detector response is taken into account
            data = calc_nll_llr_dtf_delta(filtered) if dtf == "delta" else (calc_nll_llr_dtf_dbgauss(filtered) if dtf == "dbgauss" else None)
            if data is None:
                raise Exception("Could not fetch likelihoods from event data. Check the transfer function")
            
            nll_zhh_path = os.path.join(dst_dir, "plots", "nll_zhh.png")
            nll_zzh_path = os.path.join(dst_dir, "plots", "nll_zzh.png")
            
            plot_nll(data, name, nll_zhh_path, nll_zzh_path)
            logger.info("Saved nll_zhh at {}".format(nll_zhh_path))
            logger.info("Saved nll_zzh at {}".format(nll_zzh_path))
            
            llr_path = os.path.join(dst_dir, "plots", "llr.png")
            plot_llr(data, name, llr_path)
            logger.info("Saved llr at {}".format(llr_path))
            
        else:
            logger.info("Skipped plotting")

@click.command("me")
@click.argument("src")
@click.argument("dst")
@click.argument("name")
@click.option("--tree", default="dataTree", help="Name of TTree in ROOT file")
@click.option("--convert", is_flag=True, default=True, help="Whether or not to newly convert the ROOT file")
@click.option("--plot", is_flag=True, default=True, help="Whether or not to plot the results")
@click.option("--calc", is_flag=True, default=True, help="Whether or not to calc the likelihoods")
@click.option("--dtf", default="delta", help="Whether or not to calc the likelihoods")
def me_command(src:str, dst:str, name:str, convert:bool, plot:bool, calc:bool = True, dtf:str = "delta"):
    return me(src, dst, name, convert, plot, calc, dtf)

@click.command("me_all")
@click.argument("src")
@click.argument("dst")
@click.option("--names", default="mcparticle,reco,truejet_matchingreco,truejet_seen,truejet_true", help="Name ROOT files")
@click.option("--convert/--no-convert", default=True)
@click.option("--plot/--no-plot", default=True)
@click.option("--calc/--no-calc", default=True)
@click.option("--dtf", default="delta", help="Name of the transfer function to use")
def me_all_command(src:str, dst:str, names:str, convert:bool, plot:bool, calc:bool, dtf:str):
    """Does the me analysis for all sub-files specified by names inside the src directory and outputs to dst/name"""
    
    for name in names.split(","):
        me(os.path.join(src, "compare_" + name + ".root"), os.path.join(dst, name), name, convert, plot, calc, dtf)