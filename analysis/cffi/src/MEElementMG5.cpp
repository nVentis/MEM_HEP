// Accessing the C++ subprocess classes through the C ABI 
// Used for both ZHH and ZZH

#include "MEElementMG5.hpp"

double MEElementMG5::single(std::vector<double*> momenta) {
    // Insert initial kinematics
    _cppp->setMomenta(momenta);
    _cppp->sigmaKin();
    return _cppp->getMatrixElements()[0]; // here: one processor per SubProcesses-folder
}

// Expects 6 dimensional phase-space (for final states), i.e. double momenta[24]
double* MEElementMG5::multi(double momenta[], int n_elements, double buffer[]) {
    // Initial states (e+,e-) are assumed to follow default format in 500 GeV simulations, i.e. E,px,py,pz=2.5e+2,1.75e+0,0,+-2.5e+2
    // Total 8*4 phase space dimensions

    // MadGraph uses E,px,py,pz convention per particle

    std::vector<double*> p(1, new double[4]);
    p[0][0] = constants.sqrt_s/2;
    p[0][1] = 0.;
    p[0][2] = 0.;
    p[0][3] = constants.sqrt_s/2;

    p.push_back(new double[4]);
    p[1][0] = constants.sqrt_s/2;
    p[1][1] = 0.;
    p[1][2] = 0.;
    p[1][3] = -constants.sqrt_s/2;

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
double MEElementMG5::rambo()
{
    double energy = constants.sqrt_s;
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