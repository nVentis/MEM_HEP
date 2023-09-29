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

# Call python me.cffi-builder.py, optionally with -vvv, -vv or -v
if __name__ == '__main__':    
    options = {
        # narrow width approximation
        'NWA': (True if (len(sys.argv) > 1 and True in [argv.lower().endswith("-nwa") for argv in sys.argv]) else False),
        
         # separate transfer functions for sig/bkg
        'SEPTF': (True if (len(sys.argv) > 1 and True in [argv.lower().endswith("-septf") for argv in sys.argv]) else False),
        
        'DEBUG_VVV': (True if (len(sys.argv) > 1 and True in [argv.lower().endswith("-vvv") for argv in sys.argv]) else False),
        'DEBUG_VV': (True if (len(sys.argv) > 1 and True in [argv.lower().endswith("-vv") for argv in sys.argv]) else False),
        'DEBUG_V': (True if (len(sys.argv) > 1 and True in [argv.lower().endswith("-v") for argv in sys.argv]) else False)
    }
    
    extra_compile_args = ['-fPIC']
    
    if options["DEBUG_VVV"]:
        print("[Compiling with DEBUG_VVV DEBUG_VV DEBUG_V macro]")
        extra_compile_args.append('-DDEBUG_VVV')
        extra_compile_args.append('-DDEBUG_VV')
        extra_compile_args.append('-DDEBUG_V')
    elif options["DEBUG_VV"]:
        print("[Compiling with DEBUG_VV DEBUG_V macro]")
        extra_compile_args.append('-DDEBUG_VV')
        extra_compile_args.append('-DDEBUG_V')
    elif options["DEBUG_V"]:
        print("[Compiling with DEBUG_V macro]")
        extra_compile_args.append('-DDEBUG_V')
        
    if options["NWA"]:
        print("[Compiling with NWA macro]")
        extra_compile_args.append('-DNWA')
        
    if options["SEPTF"]:
        print("[Compiling with SEPTF macro]")
        
    f = open("./compiled_with.py", "w")
    f.write(
'''lib_options = {
''' + ("\n".join([ f"   '{option}': {(str('True') if options[option] else str('False'))}," for option in options ])) + '''
}''')
    f.close()
    
    source = f'''double calc_rambo(const char param_card[], double energy, int helicity_selection[], int selected_helicities);
                double* calc(const char param_card[], double energy, int helicity_selection[], int selected_helicities, double momenta[], int n_elements);
                double* calc_mc_batch(const char param_card[], double energy, int helicity_selection[], int selected_helicities, double reco_kin[], double int_variables[], int n_elements, int calc_mc_batch);
                int calc_kinematics_from_int(const char param_card[], double evt_constants[], int helicity_selection[], int selected_helicities, {'double mH2,' if not options["NWA"] else '' } double Thb1, double Phb1, double Rhb1, double Thb1b, double Phb1b, double Rhb2, double Thb2);'''
    
    for process in processes:
        process_path = processes[process]
        extra_compile_args_proc = []
        if options["SEPTF"]:
            extra_compile_args_proc.append("-D" + ("SIG" if process == "zhh" else "BKG") + "HYP")
        
        extra_compile_args_run = extra_compile_args + extra_compile_args_proc
        
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
            extra_compile_args=extra_compile_args_run
        )

        # create C code for module and compile it
        ffibuilder.compile(verbose=True)
        