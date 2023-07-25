import sys
import os
import shutil
# Expects two arguments:
# 1. comparison_root file
# 2. comparison_out directory: Where to store the converted data and plots
from ..analysis.convert_directory import convert_file
from ..analysis.import_data import import_data, filter_data
from ..analysis.calc import calc_nll

if __name__ == '__main__':

    src_file = sys.argv[1]
    dst_dir  = sys.argv[2]
    prod_plt_mode = sys.argv[3] if len(sys.argv) > 3 else 0

    if prod_plt_mode:
        from ..analysis.plot_root import plot_hist
    else:
        from ..analysis.plot_matplotlib import plot_hist

    filename = os.path.basename(src_file)

    cnv_file = os.path.join(dst_dir, os.path.splitext(filename)[0] + ".npy")
    for file in [cnv_file]:
        if os.path.isfile(file):
            os.remove(file)

    plots_dir = os.path.join(dst_dir, "plots")
    for dir in [plots_dir]:
        if os.path.isdir(dir):
            shutil.rmtree(dir, ignore_errors=True)
            
        os.makedirs(dir)

    convert_file(src_file, "dataTree", cnv_file)