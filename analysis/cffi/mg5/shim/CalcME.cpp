// Accessing the C++ subprocess classes through the C ABI 
// Used for both ZHH and ZZH

#include "CalcME.hpp"

//----------- Implementation of class calc_me -----------//

calc_me::calc_me(std::string param_card, double energy)
{
    _cppp = new CPPProcess();
    constants.sqrt_s = energy;

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

double calc_me::calc_jac(
    double Rhb1, double Thb1, double Phb1,
    double Rhb1b, double Thb1b, double Phb1b,
    double Rhb2, double Thb2, double Phb2,
    double Rhb2b, double Thb2b, double Phb2b
){

    return std::abs(
        (1/std::sqrt(mb_pow2 + std::pow(Rhb2b, 2.) )) * (
            Rhb2*(std::pow(Rhb2b, 3.)*std::sin(Phb2 - Phb2b)* 2*(
                + std::sqrt(mb_pow2 + std::pow(Rhb1, 2.) )*Rhb1b/std::sqrt(mb_pow2 + std::pow(Rhb1b, 2.))
                - Rhb1*(
                      std::cos(Thb1)*std::cos(Thb1b)
                    + std::cos(Phb1 - Phb1b)*std::sin(Thb1)*std::sin(Thb1b)
                )
            )*std::sin(Thb2)*std::pow(std::sin(Thb2b), 2.))
        )
    );
}

int calc_me::calc_kinematics_from_int(
    double mH2, double Thb1, double Phb1, double Rhb1,
    double Thb1b, double Phb1b, double Rhb2, double Thb2
) {
    double arg_b1E_sqrt = mb_pow2 + std::pow(Rhb1, 2.);
    if (arg_b1E_sqrt < 0)
        return err_map[0];
    
    double b1E = std::sqrt(arg_b1E_sqrt);
    
    vec3 b1e = vec3_sph(Thb1, Phb1, 1);
    
    // Calculate Rhb1b from mH2, Rhb1 and angles
    vec3 b1be = vec3_sph(Thb1b, Phb1b, 1);
    double dot = vec3_dot(b1e, b1be);
    
    double dmH2 = (mH2 - 2*mb_pow2)/2;
    double d = (dmH2*Rhb1*dot)/(std::pow(Rhb1*dot, 2.) - std::pow(b1E, 2.) );
    
    double arg_Rhb1b_sqrt = (std::pow(b1E, 2.)*mb_pow2 - std::pow(dmH2, 2.))/(std::pow(Rhb1*dot, 2.) - std::pow(b1E, 2.)) + std::pow(d, 2.);
    if (arg_Rhb1b_sqrt < 0)
        return err_map[1];

    double Rhb1b = -d + std::sqrt(arg_Rhb1b_sqrt);
    if (Rhb1b < 0)
        return err_map[2];
    
    double arg_b1bE_sqrt = mb_pow2 + std::pow(Rhb1b, 2.);
    if (arg_b1bE_sqrt < 0)
        return err_map[3];
    
    double b1bE = std::sqrt(arg_b1bE_sqrt);
    
    std::array<double, 4> pB2 {
        constants.system_E -(b1E+b1bE),
        constants.system_px -Rhb1*b1e[0] -Rhb1b*b1be[0],
        constants.system_py -Rhb1*b1e[1] -Rhb1b*b1be[1],
        constants.system_pz -Rhb1*b1e[2] -Rhb1b*b1be[2],
    };

    // Calculate Rhb2b
    double arg_b2E_sqrt = mb_pow2 + std::pow(Rhb2, 2.);
    if (arg_b2E_sqrt < 0)
        return err_map[4];
    double b2E = sqrt(arg_b2E_sqrt);
    
    double b2bE = pB2[0] - b2E;
    if (b2bE < 0)
        return err_map[5];
    
    double arg_Rhb2b_sqrt = std::pow(b2bE, 2.) - mb_pow2;
    if (arg_Rhb2b_sqrt < 0)
        return err_map[6];
    double Rhb2b = std::sqrt(arg_Rhb2b_sqrt);
    
    // Calculate remaining variables, i.e. Thb2b, Phb2b and Phb2
    double arg_Thb2b_acos = 1/Rhb2b*(pB2[3] - Rhb2*std::cos(Thb2));
    if (std::abs(arg_Thb2b_acos) > 1)
        return err_map[7];
    double Thb2b = acos(arg_Thb2b_acos);
    
    double a_pow2 = std::pow(pB2[1], 2.);
    double b_pow2 = std::pow(pB2[2], 2.);
    double c = (-std::pow(Rhb2*sin(Thb2), 2.) + std::pow(Rhb2b*sin(Thb2b), 2.) + a_pow2 + b_pow2)/(2*Rhb2b*std::sin(Thb2b));
    
    double arg_Phb2b_atan_sqrt = a_pow2 + b_pow2 - std::pow(c, 2.);
    if (arg_Phb2b_atan_sqrt < 0)
        if (-arg_Phb2b_atan_sqrt < constants.epsilon) {
            arg_Phb2b_atan_sqrt = 0;
        } else {
            return err_map[8];
        }

    // Two solutions to calculate Phb2b with atan; try out both, use the one that is closest to target   
    double Phb2b_1 = 2*std::atan((pB2[2] + sqrt(arg_Phb2b_atan_sqrt))/(pB2[1]+c));
    double Phb2b_2 = 2*std::atan((pB2[2] - sqrt(arg_Phb2b_atan_sqrt))/(pB2[1]+c));
    
    double arg_Phb2_1_acos = (pB2[1] - Rhb2b*std::sin(Thb2b)*std::cos(Phb2b_1))/(Rhb2*sin(Thb2));
    double arg_Phb2_2_acos = (pB2[1] - Rhb2b*std::sin(Thb2b)*std::cos(Phb2b_2))/(Rhb2*sin(Thb2));
    
    if (std::abs(arg_Phb2_1_acos) && std::abs(arg_Phb2_2_acos) > 1)
        return -err_map[9];
    
    // Check which solution is better    
    double Phb2 = 0;
    double Phb2b = 0;
    if (std::abs(arg_Phb2_1_acos) > 1) {
        Phb2b = Phb2b_2;
        Phb2 = std::acos(arg_Phb2_2_acos);
        
        if ((pB2[2] - Rhb2b*std::sin(Thb2b)*std::sin(Phb2b))/(std::sin(Thb2)) < 0)
            Phb2 = -Phb2;
    } else if (std::abs(arg_Phb2_2_acos) > 1) {
        Phb2b = Phb2b_1;
        Phb2 = std::acos(arg_Phb2_1_acos);
        
        if ((pB2[2] - Rhb2b*std::sin(Thb2b)*std::sin(Phb2b))/(std::sin(Thb2)) < 0)
            Phb2 = -Phb2;
    } else {
        double Phb2_1 = acos(arg_Phb2_1_acos);
        double Phb2_2 = acos(arg_Phb2_2_acos);
        
        if ((pB2[2] - Rhb2b*sin(Thb2b)*sin(Phb2b_1))/sin(Thb2) < 0)
            Phb2_1 = -Phb2_1;
            
        if ((pB2[2] - Rhb2b*sin(Thb2b)*sin(Phb2b_2))/sin(Thb2) < 0)
            Phb2_2 = -Phb2_2;
        
        vec3 b2_1p = vec3_sph(Thb2, Phb2_1, Rhb2);
        vec3 b2_2p = vec3_sph(Thb2, Phb2_2, Rhb2);

        vec3 b2b_1p = vec3_sph(Thb2b, Phb2b_1, Rhb2b);
        vec3 b2b_2p = vec3_sph(Thb2b, Phb2b_2, Rhb2b);
        
        vec3 B1check_1 {
            pB2[1] - b2_1p[0] - b2b_1p[0],
            pB2[2] - b2_1p[1] - b2b_1p[1],
            pB2[3] - b2_1p[2] - b2b_1p[2]
        };

        vec3 B1check_2 {
            pB2[1] - b2_2p[0] - b2b_2p[0],
            pB2[2] - b2_2p[1] - b2b_2p[1],
            pB2[3] - b2_2p[2] - b2b_2p[2]
        };
        
        if (vec3_dot(B1check_1, B1check_1) < vec3_dot(B1check_2, B1check_2)) {
            Phb2 = Phb2_1;
            Phb2b = Phb2b_1;
        } else {
            Phb2 = Phb2_2;
            Phb2b = Phb2b_2;
        }
    }

    // Finish
    vec3 b2p = vec3_sph(Thb2, Phb2, Rhb2);
    vec3 b2bp = vec3_sph(Thb2b, Phb2b, Rhb2b);
    
    // First spherical coordinates
    kin[0] = Rhb1;
    kin[1] = Thb1;
    kin[2] = Phb1;
    
    kin[3] = Rhb1b;
    kin[4] = Thb1b;
    kin[5] = Phb1b;

    kin[6] = Rhb2;
    kin[7] = Thb2;
    kin[8] = Phb2;

    kin[9] = Rhb2b;
    kin[10] = Thb2b;
    kin[11] = Phb2b;

    // Second, four vectors
    kin[12] = b1E;
    kin[13] = Rhb1*b1e[0];
    kin[14] = Rhb1*b1e[1];
    kin[15] = Rhb1*b1e[2];

    kin[16] = b1bE;
    kin[17] = Rhb1*b1be[0];
    kin[18] = Rhb1*b1be[1];
    kin[19] = Rhb1*b1be[2];

    kin[20] = b2E;
    kin[21] = b2p[0];
    kin[22] = b2p[1];
    kin[23] = b2p[2];

    kin[24] = b2bE;
    kin[25] = b2bp[0];
    kin[26] = b2bp[1];
    kin[27] = b2bp[2];

    // Check validity
    std::array<double, 4> check {
        constants.system_E -kin[12] -kin[16] -kin[20] -kin[24],
        constants.system_px -kin[13] -kin[17] -kin[21] -kin[25],
        constants.system_py -kin[14] -kin[18] -kin[22] -kin[26],
        constants.system_pz -kin[15] -kin[19] -kin[23] -kin[27]
    };
    
    double dE = check[0];
    double dp = sqrt(std::pow(check[1], 2) + std::pow(check[2], 2) + std::pow(check[3], 2));
    
    if (dE > constants.dEmax || dp > constants.dpmax)
        return err_map[10];
    else
        return 1;
}

void calc_me::mem_init(double evt_constants[]) {
    constants.mb = evt_constants[0];
    constants.epsilon = evt_constants[1];
    constants.dEmax = evt_constants[2];
    constants.dpmax = evt_constants[3];
    constants.sqrt_s = evt_constants[4];

    constants.system_E = evt_constants[5]; // i.e. sqrt_s - mu1E -mu2E
    constants.system_px = evt_constants[6]; // same for momenta
    constants.system_py = evt_constants[7];
    constants.system_pz = evt_constants[8];
}

void calc_me::mc_batch(double reco_kin[], double int_variables[], int n_elements, double buffer[]) {
    // reco_kin: 24-element array (E,px,py,pz) for (mu-,mu+,b1,b1b,b2,b2b)
    // int_variables: array of length n_elements*8 
    // n_elements: must match to int_variables and buffer
    // buffer: array of length n_elements
    
    int kin_result = 0.;
    double me_element = 0.;
    double transfer = 0.;
    double jacobian = 0.;

    int i,j,k;

    // 
    std::array<std::array<double, 2>, 4> reco_angles = get_reco_angles(reco_kin);

    constants.system_E = constants.sqrt_s - reco_kin[0] - reco_kin[4];
    constants.system_px =               0 - reco_kin[1] - reco_kin[5];
    constants.system_py =               0 - reco_kin[2] - reco_kin[6];
    constants.system_pz =               0 - reco_kin[3] - reco_kin[7];

    // Initialize full_kin compliant to MG5
    // Set initial e-e+ momenta
    std::vector<double*> full_kin(1, new double[4]);
    full_kin[0][0] = constants.sqrt_s/2;
    full_kin[0][1] = 0.;
    full_kin[0][2] = 0.;
    full_kin[0][3] = constants.sqrt_s/2;

    full_kin.push_back(new double[4]);
    full_kin[1][0] = constants.sqrt_s/2;
    full_kin[1][1] = 0.;
    full_kin[1][2] = 0.;
    full_kin[1][3] = -constants.sqrt_s/2;

    // Set muon momenta
    full_kin.push_back(new double[4]);
    full_kin.push_back(new double[4]);

    for (i = 0; i < 2; i++) {
        for (j = 0; j < 4; j++) {
            full_kin[2+i][j] = reco_kin[(i*4) + j];
        }
    }

    // Prepare vectors for quar quantities
    full_kin.push_back(new double[4]);
    full_kin.push_back(new double[4]);
    full_kin.push_back(new double[4]);
    full_kin.push_back(new double[4]);

    for (i = 0; i < n_elements; i++) {
        if(calc_kinematics_from_int(
            int_variables[i*8],
            int_variables[i*8+1],
            int_variables[i*8+2],
            int_variables[i*8+3],
            int_variables[i*8+4],
            int_variables[i*8+5],
            int_variables[i*8+6],
            int_variables[i*8+7]
        ) < 0){
            buffer[i] = 0.;
        } else {
            // Adjust current parton four-momenta, evaluate matrix element
            for (j = 0; j < 4; j++) {
                for (k = 0; k < 4; k++) {
                    full_kin[2+j][k] = kin[12+j];
                }
            }

            me_element = calc_me_single_full_kin(full_kin);
            
            // Energy and angle transfer
            transfer = 1.;
            for (j = 0; j < 4; j++) {
                transfer = transfer*calc_tf_E(reco_kin[8+j*4], kin[12 + j*4]);
                transfer = transfer*calc_tf_Th(reco_angles[j][0], kin[1 + j*3]);
                transfer = transfer*calc_tf_Ph(reco_angles[j][1], kin[2 + j*3]);
            }

            // Calculate jacobian
            jacobian = calc_jac(
                kin[0], kin[1], kin[2],
                kin[3], kin[4], kin[5],
                kin[6], kin[7], kin[8],
                kin[9], kin[10], kin[11]
            );
                
            buffer[i] = me_element*transfer*jacobian;
        }
    }

    // For every new, a delete...
    for (double* pointer : full_kin) {
        delete pointer;
    }
    full_kin.clear();
}

double calc_me::calc_tf_E(double a, double b) {
    // a: Measured, b:True
    return 1/(M_PI*tf_E_args[1]*(1 + std::pow( ((a-b)-tf_E_args[0])/tf_E_args[1], 2.) ));
}

double calc_me::calc_tf_Th(double a, double b) {
    // a:Measured, b:True
    return 1/(M_PI*tf_Th_args[1]*(1 + std::pow( ((a-b)-tf_Th_args[0])/tf_Th_args[1], 2.) ));
}

double calc_me::calc_tf_Ph(double a, double b) {
    // a:Measured, b:True
    return 1/(M_PI*tf_Ph_args[1]*(1 + std::pow( ((a-b)-tf_Ph_args[0])/tf_Ph_args[1], 2.) ));
}

double calc_me::calc_me_single_full_kin(std::vector<double*> momenta) {
    // Insert initial kinematics
    _cppp->setMomenta(momenta);
    _cppp->sigmaKin();
    return _cppp->getMatrixElements()[0]; // here: one processor per SubProcesses-folder
}

// Expects 6 dimensional phase-space (for final states), i.e. double momenta[24]
double* calc_me::calc_me_multi(double momenta[], int n_elements, double buffer[]) {
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

double* calc_me_mc_batch(pStat self, double reco_kin[], double int_variables[], int n_elements)
{
    double *buffer;
    buffer = (double*)malloc(n_elements * sizeof(double));

    auto p = reinterpret_cast<calc_me*>(self);
    p->mc_batch(reco_kin, int_variables, n_elements, buffer);

    return buffer;
}

double calc_me_rambo(pStat self)
{
    auto p = reinterpret_cast<calc_me*>(self);
    return p->calc_me_rambo();
}
