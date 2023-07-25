import sys
import os
import shutil
import logging
logger = logging.getLogger('mem_hep')

# Expects two arguments:
# 1. comparison_root file
# 2. comparison_out directory: Where to store the converted data and plots
# 3. optional: production plot mode; if true, uses ROOT, otherwise uses matplotlib

from analysis.convert_directory import convert_file
from analysis.import_data import import_data, filter_data, combine_columns, split_true_zhh_zzh
from analysis.calc import calc_nll_llr_dtf_delta
from analysis.plot_routines import plot_summary, plot_nll, plot_llr

def mem_plot(src_file, dst_dir, name, convert, plot):
    logger.info("Analyzing {}".format(name))
    
    root_file = os.path.basename(src_file)
    cnv_file = os.path.join(dst_dir, os.path.splitext(root_file)[0] + ".npy")
    
    if convert:
        for file in [cnv_file]:
            if os.path.isfile(file):
                logger.info("Deleting file {}".format(file))
                os.remove(file)
                
        if not os.path.isdir(dst_dir):
            os.makedirs(dst_dir)
            logger.info("Created directory tree {}".format(dst_dir))
                
        convert_file(src_file, "dataTree", cnv_file)
        logger.info("Converted file {} to {}".format(src_file, cnv_file))
    else:
        logger.info("Skipped conversion")

    if plot:
        plots_dir = os.path.join(dst_dir, "plots")
        for dir in [plots_dir]:
            if os.path.isdir(dir):
                shutil.rmtree(dir, ignore_errors=True)
                logger.info("Clearing directory-tree {}".format(dir))
                
            os.makedirs(dir)
            
        raw = import_data(cnv_file)

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
        data = calc_nll_llr_dtf_delta(filtered)
        
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
                
        
        
        