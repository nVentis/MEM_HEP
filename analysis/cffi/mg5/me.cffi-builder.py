# Compile in py311 environment
# Add '-fPIC' to make files in mg5 necessary for getting rambo working!!!

from cffi import FFI
from os import path as osp

# Compiles all sub-processes with the same shim files
processes = {
    'zhh': 'P1_Sigma_sm_emep_mummupbbxbbx',
    'zzh': 'P2_Sigma_sm_emep_mummupbbxbbx'
}

if __name__ == '__main__':
    for process in processes:
        process_path = processes[process]
        
        ffibuilder = FFI()

        # specify functions, etc to be made available to Python
        ffibuilder.cdef('''double calc_rambo(const char param_card[], double energy, int helicity_selection[], int selected_helicities);
                        double* calc(const char param_card[], double energy, int helicity_selection[], int selected_helicities, double momenta[], int n_elements);
                        void free(void *ptr);''')

        # specify code needed to build the Python module
        ffibuilder.set_source(
            module_name=f'CalcME{process.upper()}',
            source='''double calc_rambo(const char param_card[], double energy, int helicity_selection[], int selected_helicities);
                    double* calc(const char param_card[], double energy, int helicity_selection[], int selected_helicities, double momenta[], int n_elements);''',
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
            extra_compile_args=['-fPIC']
        )

        # create C code for module and compile it
        ffibuilder.compile(verbose=True)
        