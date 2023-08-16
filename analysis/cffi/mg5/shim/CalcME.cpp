// Accessing the C++ subprocess classes through the C ABI 
// Used for both ZHH and ZZH

#include "CalcME.hpp"

//----------- Implementation of class calc_me -----------//

calc_me::calc_me(std::string param_card, double energy)
{
    _cppp = new CPPProcess();
    E = energy;

    //std::cout << _cppp->name() << std::endl;

    _cppp->initProc(param_card);

    //std::cout << " [TRACE] CPPP created Ok." << std::endl;
}

calc_me::~calc_me()
{
    //std::cout << " [TRACE] CPPP deleted OK" << std::endl;
}

void calc_me::set_helicity(int particle, int helicity)
{
    //std::cout << " [TRACE] set_helicity" << particle << "->" << helicity << std::endl;
    _cppp->setHelicity(particle, helicity);
}

// Expects 6 dimensional phase-space (for final states), i.e. double momenta[24]
double* calc_me::calc_me_multi(double momenta[], int n_elements, double buffer[]) {
    // Initial states (e+,e-) are assumed to follow default format in 500 GeV simulations, i.e. E,px,py,pz=2.5e+2,1.75e+0,0,+-2.5e+2
    // Total 8*4 phase space dimensions

    // MadGraph uses E,px,py,pz convention per particle

    std::vector<double*> p(1, new double[4]);
    p[0][0] = E/2;
    p[0][1] = 0.;
    p[0][2] = 0.;
    p[0][3] = E/2;

    p.push_back(new double[4]);
    p[1][0] = E/2;
    p[1][1] = 0.;
    p[1][2] = 0.;
    p[1][3] = -E/2;

    p.push_back(new double[4]);
    p.push_back(new double[4]);
    p.push_back(new double[4]);
    p.push_back(new double[4]);
    p.push_back(new double[4]);
    p.push_back(new double[4]);
    
    int k = 0;
    for (int j = 0; j < n_elements; j++) {
        // Update final state momenta
        for(int i=0; i < 6; i++){
            p[2+i][0] = momenta[k*24 + 4*i];
            p[2+i][1] = momenta[k*24 + 4*i+1];
            p[2+i][2] = momenta[k*24 + 4*i+2];
            p[2+i][3] = momenta[k*24 + 4*i+3];
        }

        _cppp->setMomenta(p);
        _cppp->sigmaKin();
        buffer[j] = _cppp->getMatrixElements()[0];

        k++;
    }

    for (double* pointer : p) {
        delete pointer;
    }
    p.clear();

    return buffer;
};

// Px,Py,Pz,E
double calc_me::calc_me_rambo()
{
    double energy = 500;
    double weight;
    
    // Get phase space point
    std::vector<double*> p = get_momenta(_cppp->ninitial, energy,
                    _cppp->getMasses(), weight);

    // Set momenta for this event
    _cppp->setMomenta(p);

    // Evaluate matrix element
    _cppp->sigmaKin();

    const double* matrix_elements = _cppp->getMatrixElements();

    return matrix_elements[0];
}

//---------- C-Interface for class calc_me ---------------------//

pStat calc_me_new(const char param_card[], double energy)
{
    return new (std::nothrow) calc_me(param_card, energy);
}

void calc_me_set_helicity(pStat self, int particle, int helicity)
{
    auto p = reinterpret_cast<calc_me*>(self);
    return p->set_helicity(particle, helicity);
}

void calc_me_del(pStat self)
{
    delete reinterpret_cast<calc_me*>(self);
}

double* calc_me_multi(pStat self, double momenta[], int n_elements)
{
    double *buffer;
    buffer = (double*)malloc(n_elements * sizeof(double));

    auto p = reinterpret_cast<calc_me*>(self);
    return p->calc_me_multi(momenta, n_elements, buffer);
}

double calc_me_rambo(pStat self)
{
    auto p = reinterpret_cast<calc_me*>(self);
    return p->calc_me_rambo();
}
