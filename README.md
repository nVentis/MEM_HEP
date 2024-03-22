# MEM_HEP

Collection of work for my Master's Thesis.

Includes 

1. code for (the namegiving) Matrix Element Method for $ZHH$ and $ZZH$ to $l\bar{l}(Z) + b\bar{b} + b\bar{b}$ including the phase-space parameterization for MC integration
    - in VEGAS
    - with a neural importance sampling strategy using rational quadratic coupling transforms by Durkan et al. [1] implemented in the normflows package by Stimper et al. [2]. The integration procedure is based on  i-flow [3].
2. the `JetConv` Marlin processor for jet clustering with GNNs and Spectral Clustering. see [GraphJet](https://gitlab.desy.de/bryan.bliewert/graphjet)

# Setup

Using [mamba and miniforge](https://github.com/conda-forge/miniforge) is recommended. Create a virtual environments

    mamba create -n py311 python=3.11 && mamba activate py311

Requires

    mamba install pytorch numpy matplotlib seaborn pandas click
    pip install vegas uproot

Alternatively, consider `conda create --name py311 --file requirements.txt` (however, this will install many modules more than what is strictly required).

To build the integrand, check `analysis/cffi` and 

For neural importane sampling, install the forked normflows package

    git clone https://github.com/nVentis/normalizing-flows.git
    pip install -e normalizing-flows

Using and converting LCIO files requires pyLCIO to be available in the conda environment. See [here](https://github.com/iLCSoft/LCIO) for build instructions.

# How-To

For the integration, check the CLI command `python cli mem_integrate` and `mem_integrate` in `analysis/mem.py`.

# Citations

[1] "Neural Spline Flows", by Durkan et al. (2019). Proceedings of the 33rd International Conference on Neural Information Processing Systems (NeurIPS 2019) [arXiv:1906.04032](https://arxiv.org/abs/1906.04032)

[2] "normflows: A PyTorch Package for Normalizing Flows", by Stimper et al. (2023). Journal of Open Source Software, 8(86), 5361, [DOI](https://doi.org/10.21105/joss.05361) [arXiv:2302.12014](https://arxiv.org/abs/2302.12014).

[3] "i-flow: High-dimensional Integration and Sampling with Normalizing Flows",
by Christina Gao, Joshua Isaacson, Claudius Krause (2020).
Mach. Learn.: Sci. Technol. [DOI](https://doi.org/10.1088/2632-2153/abab62) 
[arXiv:2001.05486](https://arxiv.org/abs/2001.05486).