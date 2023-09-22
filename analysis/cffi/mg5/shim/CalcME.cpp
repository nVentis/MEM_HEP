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

double calc_me::calc_jac(double params[]){
    /*
    Rhb1 = params[0]
    Thb1 = params[1]
    Phb1 = params[2]

    Rhb1b = params[3]
    Thb1b = params[4]
    Phb1b = params[5]

    Rhb2 = params[6]
    Thb2 = params[7]
    Phb2 = params[8]

    Rhb2b = params[9]
    Thb2b = params[10]
    Phb2b = params[11] */

    return std::abs(
        (1/std::sqrt(mb_pow2 + std::pow(params[9], 2.) )) * (
            params[6]*(std::pow(params[9], 3.)*std::sin(params[8] - params[11])* 2*(
                + std::sqrt(mb_pow2 + std::pow(params[0], 2.) )*params[3]/std::sqrt(mb_pow2 + std::pow(params[3], 2.))
                - params[0]*(
                      std::cos(params[1])*std::cos(params[4])
                    + std::cos(params[2] - params[5])*std::sin(params[1])*std::sin(params[4])
                )
            )*std::sin(params[7])*std::pow(std::sin(params[10]), 2.))
        )
    );
}

double calc_me::get_kinematics_from_int(double int_variables[]) {
    double mH2 = int_variables[0];
    double Thb1 = int_variables[1];
    double Phb1 = int_variables[2];
    double Rhb1 = int_variables[3];
    double Thb1b = int_variables[4];
    double Phb1b = int_variables[5];
    double Rhb2 = int_variables[6];
    double Thb2 = int_variables[7];

    double arg_b1E_sqrt = mb_pow2 + std::pow(Rhb1, 2.);
    if (arg_b1E_sqrt < 0)
        return err_map[0];
    
    double b1E = std::sqrt(arg_b1E_sqrt);
    
    vec3 b1e = spherical_vec(Thb1, Phb1, 1);
    
    // Calculate Rhb1b from mH2, Rhb1 and angles
    vec3 b1be = spherical_vec(Thb1b, Phb1b, 1);
    double dot = vec3_dot(b1e, b1be);
    
    double dmH2 = (mH2 - 2*mb_pow2)/2;
    double d = (dmH2*Rhb1*dot)/(std::pow(Rhb1*dot, 2.) - std::pow(b1E, 2.) );
    
    double arg_Rhb1b_sqrt = (std::pow(b1E, 2.)*mb_pow2 - std::pow(dmH2, 2.))/(std::pow(Rhb1*dot, 2.) - std::pow(b1E, 2.)) + std::pow(d, 2.);
    if (arg_Rhb1b_sqrt < 0)
        return err_map[1];

    double Rhb1b = -d + std::sqrt(arg_Rhb1b_sqrt);
    if (Rhb1b < 0)
        return err_map[2];
    
    double arg_b1bE1_sqrt = mb_pow2 + std::pow(Rhb1b, 2.);
    if (arg_b1bE1_sqrt < 0)
        return err_map[3];
    
    double b1bE1 = std::sqrt(arg_b1bE1_sqrt);
    
    std::array<double, 4> pB2 {
        constants["sqrt_s"] -(b1E+b1bE1),
        constants["system_px"] -b1p[0] -b1bp1[0],
        constants["system_py"] -b1p[1] -b1bp1[1],
        constants["system_pz"] -b1p[2] -b1bp1[2],
    };
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
