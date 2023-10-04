# Compile in py311 environment
# Add '-fPIC' to make files in mg5 necessary for getting rambo working!!!

from cffi import FFI
from os import path as osp

ffibuilder = FFI()

# specify functions, etc to be made available to Python
ffibuilder.cdef('''double calc_zhh_single(double momenta[]);
                double calc_zzh_single(double momenta[]);''')

# specify code needed to build the Python module
ffibuilder.set_source(
    module_name='CalculateME',
    source='''double calc_zhh_single(double momenta[]);
                double calc_zzh_single(double momenta[]);''',
    sources=[
        './CalculateME/main.c',
        './CalculateME/CalculateME.cpp'
    ],     # other sources -- file containing fcn(x, dim)
    library_dirs=[
        '/nfs/dust/ilc/user/bliewert/.mambaforge/envs/py311/lib',
        '/afs/desy.de/user/b/bliewert/public/ILCSoft/Physsim/lib64',
        
    ],
    include_dirs=[
        '/afs/desy.de/user/b/bliewert/public/ILCSoft/Physsim/include', # Physsim
        '/nfs/dust/ilc/user/bliewert/.mambaforge/envs/py311/include'
    ],
    libraries=['m', 'dl', 'Physsim'],
    extra_compile_args='-pthread -m64 -Wl,-rpath,/nfs/dust/ilc/user/bliewert/.mambaforge/envs/py311/lib -rdynamic'.split(' ')
)

if __name__ == '__main__':
    # create C code for module and compile it
    ffibuilder.compile(verbose=True)