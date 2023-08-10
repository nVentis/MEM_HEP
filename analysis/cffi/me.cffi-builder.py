from cffi import FFI
from os import path as osp

ffibuilder = FFI()
PATH = osp.dirname(__file__)

# specify functions, etc to be made available to Python
#ffibuilder.cdef('double fcn(double x[], int dim);')
ffibuilder.cdef('double calc();')

# specify code needed to build the Python module
ffibuilder.set_source(
    module_name='CalculateME',
    #source='double fcn(double x[], int dim);',
    #header_source='',
    source='double calc();',
    sources=[
        '../../source/CalculateME/main.c',
        '../../source/CalculateME/CalculateME.cpp'
    ],     # other sources -- file containing fcn(x, dim)
    library_dirs=[
        '/nfs/dust/ilc/user/bliewert/.mambaforge/envs/py311/lib',
        '/afs/desy.de/user/b/bliewert/public/ILCSoft/Physsim/lib64'
    ],
    include_dirs=[
        '/afs/desy.de/user/b/bliewert/public/ILCSoft/Physsim/include', # Physsim
        '/nfs/dust/ilc/user/bliewert/.mambaforge/envs/py311/include'
    ],
    libraries=['m', 'dl', 'Physsim'],        # may need to specify the math library (-lm)
    # 'Core', 'Imt', 'RIO', 'Net', 'Hist', 'Graf', 'Graf3d', 'Gpad', 'ROOTVecOps', 'Tree', 'TreePlayer', 'Rint', 'Postscript', 'Matrix', 'Physics', 'MathCore', 'Thread', 'MultiProc', 'ROOTDataFrame'
    # output of root-config --cflags --libs
    extra_compile_args='-pthread -m64 -Wl,-rpath,/nfs/dust/ilc/user/bliewert/.mambaforge/envs/py311/lib -rdynamic'.split(' ')
)

# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/afs/desy.de/user/b/bliewert/public/ILCSoft/Physsim/lib64

if __name__ == '__main__':
    # create C code for module and compile it
    ffibuilder.compile(verbose=True)