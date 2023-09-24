# Compile in py311 environment
# Add '-fPIC' to make files in mg5 necessary for getting rambo working!!!

from cffi import FFI
from os import path as osp
import sys

# Compiles all sub-processes with the same shim files
processes = {
    'zhh': 'P1_Sigma_sm_emep_mummupbbxbbx',
    'zzh': 'P2_Sigma_sm_emep_mummupbbxbbx'
}

# Call python me.cffi-builder.py, optionally with -vvv
if __name__ == '__main__':
    debug_vvv = True if (len(sys.argv) > 1 and "vvv" in sys.argv[1].lower()) else False
    debug_v   = True if (len(sys.argv) > 1 and "v" in sys.argv[1].lower()) else False
    
    extra_compile_args = ['-fPIC']
    
    if debug_vvv:
        print("[Compiling with DEBUG_VVV macro]")
        extra_compile_args.append('-DDEBUG_V')
        extra_compile_args.append('-DDEBUG_VVV')
    elif debug_v:
        print("[Compiling with DEBUG_V macro]")
        extra_compile_args.append('-DDEBUG_V')
    
    source = '''double calc_rambo(const char param_card[], double energy, int helicity_selection[], int selected_helicities);
                double* calc(const char param_card[], double energy, int helicity_selection[], int selected_helicities, double momenta[], int n_elements);
                double* calc_mc_batch(const char param_card[], double energy, int helicity_selection[], int selected_helicities, double reco_kin[], double int_variables[], int n_elements);
                int calc_kinematics_from_int(const char param_card[], double evt_constants[], int helicity_selection[], int selected_helicities, double mH2, double Thb1, double Phb1, double Rhb1, double Thb1b, double Phb1b, double Rhb2, double Thb2);'''
    
    for process in processes:
        process_path = processes[process]
        
        ffibuilder = FFI()

        # specify functions, etc to be made available to Python
        ffibuilder.cdef(f'''{source}
                        void free(void *ptr);''')

        # specify code needed to build the Python module
        ffibuilder.set_source(
            module_name=f'CalcME{process.upper()}',
            source=source,
            sources=[
                './shim/main.c',
                f'./mg5/SubProcesses/{process_path}/CPPProcess.cc',
                './mg5/src/Parameters_sm.cc',
                './mg5/src/read_slha.cc',
                './mg5/src/HelAmps_sm.cc',
                './mg5/src/rambo.cc',
                './shim/CalcME.cpp',
            ],
            library_dirs=[
                './mg5/lib',
            ],
            include_dirs=[
                './mg5/src', # zhh process
                f'./mg5/SubProcesses/{process_path}'
            ],
            libraries=['m', 'model_sm'],
            extra_compile_args=extra_compile_args
        )

        # create C code for module and compile it
        ffibuilder.compile(verbose=True)
        