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

# default: python me.cffi-builder.py -nwa -septf
# Call python me.cffi-builder.py, optionally with -vvv, -vv or -v
if __name__ == '__main__':
    pyEnv = 'py311'
    
    options = {
        # narrow width approximation
        'NWA': (True if (len(sys.argv) > 1 and True in [argv.lower().endswith("-nwa") for argv in sys.argv]) else False),
        
         # separate transfer functions for sig/bkg
        'SEPTF': (True if (len(sys.argv) > 1 and True in [argv.lower().endswith("-septf") for argv in sys.argv]) else False),
        
        'DEBUG_VVV': (True if (len(sys.argv) > 1 and True in [argv.lower().endswith("-vvv") for argv in sys.argv]) else False),
        'DEBUG_VV': (True if (len(sys.argv) > 1 and True in [argv.lower().endswith("-vv") for argv in sys.argv]) else False),
        'DEBUG_V': (True if (len(sys.argv) > 1 and True in [argv.lower().endswith("-v") for argv in sys.argv]) else False),
    }
    
    extra_compile_args = [
        '-fPIC',
        '-pthread',
        #'-m64',
        #'-Wl,-rpath,/nfs/dust/ilc/user/bliewert/.mambaforge/envs/py311/lib',
        #'-rdynamic'
        '-shared'
    ]
    
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
        
    # Save compile options so library can check format
    f = open("./compiled_with.py", "w")
    f.write(
'''lib_options = {
''' + ("\n".join([ f"   '{option}': {(str('True') if options[option] else str('False'))}," for option in options ])) + '''
}''')
    f.close()
    
    source = f'''double calc_rambo(const char param_card[], double energy, int helicity_selection[], int selected_helicities);
                double* calc(const char param_card[], double energy, int helicity_selection[], int selected_helicities, double momenta[], int n_elements);
                double* calc_mc_batch(const char param_card[], double energy, int helicity_selection[], int selected_helicities, double reco_kin[], double int_variables[], int n_elements, int use_transer_funcs, int me_type);
                double* calc_mc_batch_with_evt_consts(const char param_card[], double energy, int helicity_selection[], int selected_helicities, double reco_kin[], double int_variables[], int n_elements, int use_transer_funcs, int me_type, double evt_constants[]);
                int calc_kinematics_from_int(const char param_card[], double evt_constants[], int helicity_selection[], int selected_helicities, {'double mH2,' if not options["NWA"] else '' } double Thb1, double Phb1, double Rhb1, double Thb1b, double Phb1b, double Rhb2, double Thb2);'''
    
    for process in processes:
        process_path = processes[process]
        extra_compile_args_proc = ["-D" +("ZHH" if process == "zhh" else "ZZH")] # for changing between PHYSSIM ZHH/ZZH MEs
        
        if options["SEPTF"]:
            extra_compile_args_proc.append("-D" + ("SIG" if process == "zhh" else "BKG") + "HYP")
        
        extra_compile_args_run = extra_compile_args + extra_compile_args_proc
        
        ffibuilder = FFI()

        # specify functions, etc to be made available to Python
        ffibuilder.cdef(f'''{source}
                        void free(void *ptr);''')

        include_dirs = [
            './mg5/src', # zhh process
            f'./mg5/SubProcesses/{process_path}',
            
            '/afs/desy.de/user/b/bliewert/public/ILCSoft/Physsim/include', # Physsim
            f'/nfs/dust/ilc/user/bliewert/miniconda3/envs/{pyEnv}/include'
        ]

        # specify code needed to build the Python module
        ffibuilder.set_source(
            module_name=f'CalcME{process.upper()}',
            source=source,
            sources=[
                f'./mg5/SubProcesses/{process_path}/CPPProcess.cc',
                './shim/main.c',
                './mg5/src/Parameters_sm.cc',
                './mg5/src/read_slha.cc',
                './mg5/src/HelAmps_sm.cc',
                './mg5/src/rambo.cc',
                './shim/CalcME.cpp'
            ],
            library_dirs=[
                './mg5/lib',
                
                f'/nfs/dust/ilc/user/bliewert/miniconda3/envs/{pyEnv}/lib',
                '/afs/desy.de/user/b/bliewert/public/ILCSoft/Physsim/lib64'
            ],
            include_dirs=include_dirs,
            libraries=['m', 'model_sm', 'dl', 'Physsim'],
            extra_compile_args=extra_compile_args_run
        )

        # create C code for module and compile it
        ffibuilder.compile(verbose=True)
        