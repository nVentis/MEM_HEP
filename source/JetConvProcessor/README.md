# Install

Requires torchlib as well as [TorchScatter](https://github.com/rusty1s/pytorch_scatter#c-api) and [TorchSparse](https://github.com/rusty1s/pytorch_sparse/#c-api). See https://github.com/pyg-team/pytorch_geometric/tree/master/examples/cpp. Note that one needs to recursively checkout TorchSparse `git submodule update --init --recursive`

    cmake -D CMAKE_INSTALL_PREFIX:PATH=/afs/desy.de/user/b/bliewert/public/DevLocal/pytorch_scatter -DCMAKE_PREFIX_PATH="..." ..

    cmake -D CMAKE_INSTALL_PREFIX:PATH=/afs/desy.de/user/b/bliewert/public/DevLocal/pytorch_sparse -DCMAKE_PREFIX_PATH="..." ..

When using a Key4hep stack (here, `2023-05-23`), to use it's PyTorch library, make sure it can be found by CMake by adding 

    export PYTHONPATH=/cvmfs/ilc.desy.de/key4hep/releases/2023-05-23/py-torch/1.13.1/x86_64-centos7-gcc12.3.0-opt/i7cg4/lib/python3.10/site-packages:$PYTHONPATH
    export LD_LIBRARY_PATH=/cvmfs/ilc.desy.de/key4hep/releases/2023-05-23/py-torch/1.13.1/x86_64-centos7-gcc12.3.0-opt/i7cg4/lib/python3.10/site-packages/torch/lib:$LD_LIBRARY_PATH
    export CMAKE_PREFIX_PATH=/cvmfs/ilc.desy.de/key4hep/releases/2023-05-23/py-torch/1.13.1/x86_64-centos7-gcc12.3.0-opt/i7cg4/lib/python3.10/site-packages/torch/share/cmake:$CMAKE_PREFIX_PATH
    export CMAKE_PREFIX_PATH=/cvmfs/ilc.desy.de/key4hep/releases/2023-05-23/py-torch/1.13.1/x86_64-centos7-gcc12.3.0-opt/i7cg4/lib/python3.10/site-packages/torch/include/torch/csrc/api:$CMAKE_PREFIX_PATH

Same for the build dependencies

    export CMAKE_PREFIX_PATH=/afs/desy.de/user/b/bliewert/public/DevLocal/pytorch_scatter/share/cmake:$CMAKE_PREFIX_PATH
    export CMAKE_PREFIX_PATH=/afs/desy.de/user/b/bliewert/public/DevLocal/pytorch_sparse/share/cmake:$CMAKE_PREFIX_PATH

For clustering based on the affinity matrix calculated with the GNN, sklearn wrapped via pybind11 is used. To build

    git clone https://github.com/pybind/pybind11.git && cd pybind11 && mkdir build && cd build && cmake -D CMAKE_INSTALL_PREFIX:PATH=/afs/desy.de/user/b/bliewert/public/DevLocal/pybind11 -DCMAKE_PREFIX_PATH="..." .. && make && make install

TODO: Use pybind from key4hep stack

