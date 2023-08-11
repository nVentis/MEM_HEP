from cffi import FFI
from os import path as osp

ffibuilder = FFI()
PATH = osp.dirname(__file__)

# specify functions, etc to be made available to Python
#ffibuilder.cdef('double fcn(double x[], int dim);')
ffibuilder.cdef('double calc(int iters);')

# Add '-fPIC' to make files in mg5 necessary for getting rambo working!!!

# specify code needed to build the Python module
ffibuilder.set_source(
    module_name='CalcME',
    source='double calc(int iters);',
    sources=[
        './shim/zhh/main.c',
        './mg5/zhh/SubProcesses/P1_Sigma_sm_emep_mummupbbxbbx/CPPProcess.cc',
        './mg5/zhh/src/Parameters_sm.cc',
        './mg5/zhh/src/read_slha.cc',
        './mg5/zhh/src/HelAmps_sm.cc',
        './mg5/zhh/src/rambo.cc',
        './shim/zhh/CalcME.cpp',
    ],
    library_dirs=[
        './mg5/zhh/lib',
    ],
    include_dirs=[
        './mg5/zhh/src', # zhh process
        './mg5/zhh/SubProcesses/P1_Sigma_sm_emep_mummupbbxbbx'
    ],
    libraries=['m', 'model_sm'],
    extra_compile_args=['-fPIC']
)

# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/afs/desy.de/user/b/bliewert/public/ILCSoft/Physsim/lib64

if __name__ == '__main__':
    # create C code for module and compile it
    ffibuilder.compile(verbose=True)