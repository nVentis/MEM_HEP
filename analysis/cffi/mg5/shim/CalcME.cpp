// Accessing the C++ subprocess classes through the C ABI 
// Used for both ZHH and ZZH

#include "CalcME.hpp"

// MEElement
MEElementPhyssimSig::MEElementPhyssimSig(PhysConstants *consts, double pol_e, double pol_p, int z_decay_mode) {
    setConstants(consts);

    _lcme = new lcme::LCMEZHH("LCMEZHH", "ZHH", constants->mH, pol_e, pol_p);
    _lcme->SetZDecayMode(z_decay_mode); // to muon; see LCME code for internal mappings to PDG ids
    _lcme->SetNoHDecay(2);
    _lcme->SetPropagator(1);
    _lcme->SetMEType(2); // 2: ME; 1: cross-sec
};

MEElementPhyssimBkg::MEElementPhyssimBkg(PhysConstants *consts, double pol_e, double pol_p, int z1_decay_mode, int z2_decay_mode) {
    setConstants(consts);

    _lcme = new lcme::LCMEZZH("LCMEZZH", "ZZH", consts->mH, pol_e, pol_p);
    _lcme->SetZDecayMode(z1_decay_mode, z2_decay_mode);
    _lcme->SetNoZDecay(0);
    _lcme->SetNoHDecay(1);
    _lcme->SetPropagator(1);
    _lcme->SetMEType(2);
};

double MEElementPhyssimSig::single(std::vector<double*> momenta) {
    /* 
    momenta[0-1]: e- e+
    momenta[2-3]: mu- mu+
    momenta[4-5]: b bbar
    momenta[6-7]: b bbar
    */
   
    TLorentzVector _lortz[4] = {
        TLorentzVector(momenta[2][1], momenta[2][2], momenta[2][3], momenta[2][0]),
        TLorentzVector(momenta[3][1], momenta[3][2], momenta[3][3], momenta[3][0]),
        TLorentzVector(momenta[4][1] + momenta[5][1],
                       momenta[4][2] + momenta[5][2],
                       momenta[4][3] + momenta[5][3],
                       momenta[4][0] + momenta[5][0]),
        TLorentzVector(momenta[6][1] + momenta[7][1],
                       momenta[6][2] + momenta[7][2],
                       momenta[6][3] + momenta[7][3],
                       momenta[6][0] + momenta[7][0])
    };

    _lcme->SetMomentumFinal(_lortz);

    return _lcme->GetMatrixElement2();
}

double MEElementPhyssimBkg::single(std::vector<double*> momenta) {
    /* 
    momenta[0-1]: e- e+
    momenta[2-3]: mu- mu+
    momenta[4-5]: b bbar
    momenta[6-7]: b bbar
    */
   
    TLorentzVector _lortz[5] = {
        TLorentzVector(momenta[2][1], momenta[2][2], momenta[2][3], momenta[2][0]),
        TLorentzVector(momenta[3][1], momenta[3][2], momenta[3][3], momenta[3][0]),
        TLorentzVector(momenta[4][1], momenta[4][2], momenta[4][3], momenta[4][0]),
        TLorentzVector(momenta[5][1], momenta[5][2], momenta[5][3], momenta[5][0]),
        TLorentzVector(momenta[6][1] + momenta[7][1],
                       momenta[6][2] + momenta[7][2],
                       momenta[6][3] + momenta[7][3],
                       momenta[6][0] + momenta[7][0])
    };

    _lcme->SetMomentumFinal(_lortz);

    return _lcme->GetMatrixElement2();
}




//----------- Implementation of class calc_me -----------//

calc_me::calc_me(std::string param_card, double energy)
{
    _cppp = new CPPProcess();
    constants.sqrt_s = energy;

    //std::cout << _cppp->name() << std::endl;

    _cppp->initProc(param_card);

    //std::cout << " [TRACE] CPPP created Ok." << std::endl;
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

#ifndef NWA
int calc_me::calc_kinematics_from_int(double mH2, double Thb1, double Phb1, double Rhb1, double Thb1b, double Phb1b, double Rhb2, double Thb2)
#else
int calc_me::calc_kinematics_from_int(double Thb1, double Phb1, double Rhb1, double Thb1b, double Phb1b, double Rhb2, double Thb2)
#endif
{  
    /**
     * Logic: Particles b1 and b1b from H decay with Hmass^2 
     * 
    */
    double b1E = std::sqrt(mb_pow2 + std::pow(Rhb1, 2.));
    
    vec3 b1e = vec3_sph(Thb1, Phb1, 1);
    
    // Calculate Rhb1b from mH2, Rhb1 and angles
    vec3 b1be = vec3_sph(Thb1b, Phb1b, 1);
    double dot = vec3_dot(b1e, b1be);
    
    double dmH2 = (
        #ifndef NWA
            mH2
        #else
            std::pow(constants.mH, 2.)
        #endif
         - 2*mb_pow2)/2;
    double d = (dmH2*Rhb1*dot)/(std::pow(Rhb1*dot, 2.) - std::pow(b1E, 2.) );
    
    double arg_Rhb1b_sqrt = (std::pow(b1E, 2.)*mb_pow2 - std::pow(dmH2, 2.))/(std::pow(Rhb1*dot, 2.) - std::pow(b1E, 2.)) + std::pow(d, 2.);
    if (arg_Rhb1b_sqrt < 0) {
        return err_map[1];
    }

    double Rhb1b = -d + std::sqrt(arg_Rhb1b_sqrt);
    if (Rhb1b < 0) {
        return err_map[2];
    }
    
    double b1bE = std::sqrt(mb_pow2 + std::pow(Rhb1b, 2.));
    
    std::array<double, 4> pB2 {
        constants.system_E -(b1E+b1bE),
        constants.system_px -Rhb1*b1e[0] -Rhb1b*b1be[0],
        constants.system_py -Rhb1*b1e[1] -Rhb1b*b1be[1],
        constants.system_pz -Rhb1*b1e[2] -Rhb1b*b1be[2],
    };

    // Calculate Rhb2b
    double b2E = std::sqrt(mb_pow2 + std::pow(Rhb2, 2.));
    
    double b2bE = pB2[0] - b2E;
    if (b2bE < 0) {
        return err_map[5];
    }
    
    double arg_Rhb2b_sqrt = std::pow(b2bE, 2.) - mb_pow2;
    if (arg_Rhb2b_sqrt < 0) {
        return err_map[6];
    }
    double Rhb2b = std::sqrt(arg_Rhb2b_sqrt);
    
    // Calculate remaining variables, i.e. Thb2b, Phb2b and Phb2
    double arg_Thb2b_acos = 1/Rhb2b*(pB2[3] - Rhb2*std::cos(Thb2));
    if (std::abs(arg_Thb2b_acos) > 1) {
        return err_map[7];
    }
    double Thb2b = std::acos(arg_Thb2b_acos);
    
    double a_pow2 = std::pow(pB2[1], 2.);
    double b_pow2 = std::pow(pB2[2], 2.);
    double c = (-std::pow(Rhb2*sin(Thb2), 2.) + std::pow(Rhb2b*std::sin(Thb2b), 2.) + a_pow2 + b_pow2)/(2*Rhb2b*std::sin(Thb2b));
    
    double arg_Phb2b_atan_sqrt = a_pow2 + b_pow2 - std::pow(c, 2.);
    if (arg_Phb2b_atan_sqrt < 0) {
        if (-arg_Phb2b_atan_sqrt < constants.epsilon) {
            arg_Phb2b_atan_sqrt = 0;
        } else {
            return err_map[8];
        }
    }

    // Two solutions to calculate Phb2b with atan; try out both, use the one that is closest to target   
    double Phb2b_1 = 2*std::atan((pB2[2] + std::sqrt(arg_Phb2b_atan_sqrt))/(pB2[1]+c));
    double Phb2b_2 = 2*std::atan((pB2[2] - std::sqrt(arg_Phb2b_atan_sqrt))/(pB2[1]+c));
    
    double arg_Phb2_1_acos = (pB2[1] - Rhb2b*std::sin(Thb2b)*std::cos(Phb2b_1))/(Rhb2*sin(Thb2));
    double arg_Phb2_2_acos = (pB2[1] - Rhb2b*std::sin(Thb2b)*std::cos(Phb2b_2))/(Rhb2*sin(Thb2));
    
    if (std::abs(arg_Phb2_1_acos) && std::abs(arg_Phb2_2_acos) > 1) {
        return err_map[9];
    }
    
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
        double Phb2_1 = std::acos(arg_Phb2_1_acos);
        double Phb2_2 = std::acos(arg_Phb2_2_acos);
        
        if ((pB2[2] - Rhb2b*std::sin(Thb2b)*std::sin(Phb2b_1))/std::sin(Thb2) < 0) {
            Phb2_1 = -Phb2_1;
        }
            
        if ((pB2[2] - Rhb2b*std::sin(Thb2b)*std::sin(Phb2b_2))/std::sin(Thb2) < 0) {
            Phb2_2 = -Phb2_2;
        }
        
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
    kin[0] = Rhb2;
    kin[1] = Thb2;
    kin[2] = Phb2;
    
    kin[3] = Rhb2b;
    kin[4] = Thb2b;
    kin[5] = Phb2b;

    kin[6] = Rhb1;
    kin[7] = Thb1;
    kin[8] = Phb1;

    kin[9] = Rhb1b;
    kin[10] = Thb1b;
    kin[11] = Phb1b;

    // Second, four vectors
    kin[12] = b2E;
    kin[13] = b2p[0];
    kin[14] = b2p[1];
    kin[15] = b2p[2];

    kin[16] = b2bE;
    kin[17] = b2bp[0];
    kin[18] = b2bp[1];
    kin[19] = b2bp[2];

    // For both ZHH and ZZH: H -> b bbar
    kin[20] = b1E;
    kin[21] = Rhb1*b1e[0];
    kin[22] = Rhb1*b1e[1];
    kin[23] = Rhb1*b1e[2];

    kin[24] = b1bE;
    kin[25] = Rhb1b*b1be[0];
    kin[26] = Rhb1b*b1be[1];
    kin[27] = Rhb1b*b1be[2];
    
    /*
    kin[12] = b1E;
    kin[13] = Rhb1*b1e[0];
    kin[14] = Rhb1*b1e[1];
    kin[15] = Rhb1*b1e[2];

    kin[16] = b1bE;
    kin[17] = Rhb1b*b1be[0];
    kin[18] = Rhb1b*b1be[1];
    kin[19] = Rhb1b*b1be[2];

    kin[20] = b2E;
    kin[21] = b2p[0];
    kin[22] = b2p[1];
    kin[23] = b2p[2];

    kin[24] = b2bE;
    kin[25] = b2bp[0];
    kin[26] = b2bp[1];
    kin[27] = b2bp[2];
    */

    // Check validity
    std::array<double, 4> check {
        constants.system_E -kin[12] -kin[16] -kin[20] -kin[24],
        constants.system_px -kin[13] -kin[17] -kin[21] -kin[25],
        constants.system_py -kin[14] -kin[18] -kin[22] -kin[26],
        constants.system_pz -kin[15] -kin[19] -kin[23] -kin[27]
    };
    
    double dp = std::sqrt(std::pow(check[1], 2.) + std::pow(check[2], 2.) + std::pow(check[3], 2.));

    #ifdef DEBUG_VVV
    // debug print
    std::cout << "dE:" << check[0] << std::endl;
    std::cout << "dp:" << dp << std::endl;

    kin_debug_print();
    #endif

    if (check[0] > constants.dEmax || dp > constants.dpmax)
        return err_map[10];
    else
        return 1;
}

void calc_me::mem_init(double evt_constants[]) {
    #ifdef DEBUG_VVV
    std::cout << "----------------------------------------------------" << std::endl;
    std::cout << "mem_init: [" << std::endl;
    std::cout << "mb        [" << evt_constants[0] << "]" << std::endl;
    std::cout << "eps       [" << evt_constants[1] << "]" << std::endl;
    std::cout << "dEmax     [" << evt_constants[2] << "]" << std::endl;
    std::cout << "dpmax     [" << evt_constants[3] << "]" << std::endl;
    std::cout << "sqrt_s    [" << evt_constants[4] << "]" << std::endl;
    std::cout << "system_E  [" << evt_constants[5] << "]" << std::endl;
    std::cout << "system_px [" << evt_constants[6] << "]" << std::endl;
    std::cout << "system_py [" << evt_constants[7] << "]" << std::endl;
    std::cout << "system_pz [" << evt_constants[8] << "]" << std::endl;
    std::cout << "----------------------------------------------------" << std::endl;
    #endif

    constants.mb = evt_constants[0];
    constants.epsilon = evt_constants[1];
    constants.dEmax = evt_constants[2];
    constants.dpmax = evt_constants[3];
    constants.sqrt_s = evt_constants[4];

    constants.system_E = evt_constants[5]; // i.e. sqrt_s - mu1E -mu2E
    constants.system_px = evt_constants[6]; // same for momenta
    constants.system_py = evt_constants[7];
    constants.system_pz = evt_constants[8];

    mb_pow2 = std::pow(constants.mb, 2.);
    mem_init_run = true;
}

void calc_me::me_set(int me_type, bool is_signal) {
    if (me_type == 0) {
        me_element = nullptr;
    } else if (me_type == 1) {
        #ifdef ZHH
        me_element = new MEElementPhyssimSig(&constants, -1., 1., 5); // 5,11 = mu-,b
        #else
        me_element = new MEElementPhyssimBkg(&constants, -1., 1., 5, 11);
        #endif
    }
}

void calc_me::mc_batch(double reco_kin[], double int_variables[], int n_elements, double buffer[], int use_transer_funcs) {
    // reco_kin: 24-element array (E,px,py,pz) for (mu-,mu+,b1,b1b,b2,b2b)
    // int_variables: array of length n_elements*8 
    // n_elements: must match to int_variables and buffer
    // buffer: array of length n_elements
    
    int kin_result = 0.;
    double me_element_val = 0.;
    double transfer = 1.;
    double tf_E = 1.;
    double tf_Th = 1.;
    double tf_Ph = 1.;

    double jacobian = 0.;

    int i,j,k;
    int kin_res = 0;

    int tf_E_order_kin[4] {20,24,12,16};
    int tf_Th_order_kin[4] {7,10,1,4};
    int tf_Ph_order_kin[4] {8,11,2,5};

    #ifdef DEBUG_VV
    std::vector<double> B1_masses;    
    #endif

    std::array<std::array<double, 2>, 4> reco_angles = get_reco_angles(reco_kin);
    #ifdef DEBUG_VVV
        std::cout << "----------------------------------------------------" << std::endl;
        std::cout << "Reco angles: [" << std::endl;
        for (i = 0; i < 4; i++) {
            std::cout << "Particle " << (i+1) << ": Th=" << reco_angles[i][0] << " | Ph=" << reco_angles[i][1] << std::endl;
        }
        std::cout << "----------------------------------------------------" << std::endl;
    #endif

    // No need to set, if already set by mem_init
    if (!mem_init_run) {
        constants.system_E = constants.sqrt_s - reco_kin[0] - reco_kin[4];
        constants.system_px =              0. - reco_kin[1] - reco_kin[5];
        constants.system_py =              0. - reco_kin[2] - reco_kin[6];
        constants.system_pz =              0. - reco_kin[3] - reco_kin[7];
    }

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
        kin_res = calc_kinematics_from_int(
            #ifndef NWA
            int_variables[i*8],
            int_variables[i*8+1],
            int_variables[i*8+2],
            int_variables[i*8+3],
            int_variables[i*8+4],
            int_variables[i*8+5],
            int_variables[i*8+6],
            int_variables[i*8+7]
            #else
            int_variables[i*7],
            int_variables[i*7+1],
            int_variables[i*7+2],
            int_variables[i*7+3],
            int_variables[i*7+4],
            int_variables[i*7+5],
            int_variables[i*7+6]
            #endif
        );

        if(kin_res < 0){
            #ifdef DEBUG_VVV
                std::cout << "Failed with error " << kin_res << std::endl;
            #endif
            buffer[i] = 0.;
        } else {
            #ifdef DEBUG_VVV
                std::cout << "Found PS-point at [" << i << "]" << std::endl;
            #endif

            #ifdef DEBUG_VV
            B1_masses.push_back(std::sqrt(
                std::pow(kin[20] + kin[24], 2.) - (
                    std::pow(kin[21] + kin[25], 2.) +
                    std::pow(kin[22] + kin[26], 2.) +
                    std::pow(kin[23] + kin[27], 2.)
                )
            ));
            #endif

            // Adjust current parton four-momenta, evaluate matrix element
            for (j = 0; j < 4; j++) {
                for (k = 0; k < 4; k++) {
                    full_kin[4+j][k] = kin[12+(j*4)+k];
                }
            }

            #ifdef DEBUG_VVV
                int i_part = 0;
                int i_dim = 0;
                std::cout << "----------------------------------------------------" << std::endl;
                std::cout << "mc_batch full_kin [" << std::endl;
                for(i_part = 0; i_part < 8; i_part++){
                    std::cout << "Particle " << int(i_part+1) << ": ";
                    for (i_dim = 0; i_dim < 4; i_dim++) {
                        std::cout << full_kin[i_part][i_dim] << " ";
                    }
                    std::cout << std::endl;
                }
                std::cout << "]" << std::endl;
                std::cout << std::endl << "----------------------------------------------------" << std::endl;
            #endif

            if (me_element == nullptr) {
                me_element_val = calc_me_single_full_kin(full_kin);
            } else {
                // For now, only support for Physsim
                me_element_val = me_element->single(full_kin);
            }
            
            // Energy and angle transfer
            // use_transer_funcs=0 for debugging
            if (use_transer_funcs > 0) {
                transfer = 1.;
                
                #ifdef DEBUG_VVV
                std::cout << "----------------------------------------------------" << std::endl;
                std::cout << "Transfer [i] -> [Energy:E_jet:E_part][Theta:Jet:Parton][Phi:Jet:Parton] -> [Transfer]" << std::endl;
                #endif

                for (j = 0; j < 4; j++) {
                    //tf_E = calc_tf_E(reco_kin[8+j*4], kin[12 + j*4]);
                    tf_Th = calc_tf_Th(reco_angles[j][0], kin[tf_Th_order_kin[j]]);
                    tf_Ph = calc_tf_Ph(reco_angles[j][1], kin[tf_Ph_order_kin[j]]);
                    tf_E = calc_tf_E(reco_kin[8+j*4], kin[tf_E_order_kin[j]]);

                    transfer = transfer*tf_Th*tf_Ph*tf_E;
                    #ifdef DEBUG_VVV
                    std::cout << "[" << (j+1) <<  "] -> [" << tf_E << ":" << reco_kin[8+j*4] << ":" << kin[tf_E_order_kin[j]] << "][" << tf_Th << ":" << reco_angles[j][0] << ":" << kin[tf_Th_order_kin[j]] << "][" << tf_Ph << ":" << reco_angles[j][1] << ":" << kin[tf_Th_order_kin[j]] << "] -> [" << transfer << "]" << std::endl;
                    #endif
                }
                #ifdef DEBUG_VVV
                std::cout << "----------------------------------------------------" << std::endl;
                #endif
            }

            // Calculate jacobian
            jacobian = calc_jac(
                kin[0], kin[1], kin[2],
                kin[3], kin[4], kin[5],
                kin[6], kin[7], kin[8],
                kin[9], kin[10], kin[11]
            );
            
            #ifdef DEBUG_VV
                std::cout << "----------------------------------------------------" << std::endl;
                std::cout << " ME:" << me_element_val;
                std::cout << " TF:" << transfer;
                std::cout << " JAC:" << jacobian << std::endl;
                std::cout << "----------------------------------------------------" << std::endl;
            #endif

            buffer[i] = me_element_val*transfer*jacobian;
        }
    }

    #ifdef DEBUG_VV
    std::cout << "----------------------------------------------------" << std::endl;
    std::cout << "B1_masses = [" << std::endl;
    for (auto x: B1_masses) {
        std::cout << x << ", ";
    }
    std::cout << std::endl << "]" << std::endl;
    std::cout << "----------------------------------------------------" << std::endl;
    #endif

    // For every new, a delete...
    for (double* pointer : full_kin) {
        delete pointer;
    }
    full_kin.clear();
}

#ifdef PARTDBGAUSS
    #ifdef DEBUG_VV
        std::cout << "Using DBGauss DTF for Partons" << std::endl;
    #endif

    double calc_me::calc_tf_E(double a, double b) {
        // a: Measured, b:True
        return 1/(std::sqrt(2*M_PI) * (tf_E_args[1] + tf_E_args[4]*std::pow(tf_E_args[2], 2.) )) * (
                                           std::exp(-std::pow((a-b)-tf_E_args[0], 2.) / (2*std::pow(tf_E_args[1], 2.)))
            + std::pow(tf_E_args[2], 2.) * std::exp(-std::pow((a-b)-tf_E_args[3], 2.) / (2*std::pow(tf_E_args[4], 2.)))
        );
    }
#else
    #ifdef DEBUG_VV
        std::cout << "Using DBGauss DTF for Partons" << std::endl;
    #endif

    double calc_me::calc_tf_E(double a, double b) {
        // a: Measured, b:True
        return 1/(M_PI*tf_E_args[1]*(1 + std::pow( ((a-b)-tf_E_args[0])/tf_E_args[1], 2.) ));
    }
#endif


double calc_me::calc_tf_Th(double a, double b) {
    // a:Measured, b:True
    return 1/(M_PI*tf_Th_args[1]*(1 + std::pow( ((a-b)-tf_Th_args[0])/tf_Th_args[1], 2.) ));
}

double calc_me::calc_tf_Ph(double a, double b) {
    // a:Measured, b:True
    return 1/(M_PI*tf_Ph_args[1]*(1 + std::pow( ((a-b)-tf_Ph_args[0])/tf_Ph_args[1], 2.) ));
}

void calc_me::kin_debug_print(){
    int t = 0;
    int s = 1;
    int m = 0;
    int n = 0;
    
    std::cout << "----------------------------------------------------" << std::endl;

    std::cout << "Spherical coordinates:" << std::endl << "Particle 1: ";
    for(auto x: kin){
        std::cout << x << " ";

        t++;
        n++;

        if (t == 27)
            break;

        if (n == 12) {
            m = 1;
            s = 1;
            t = 0;
            std::cout << std::endl << std::endl << "Four-vectors:" << std::endl << "Particle 1: ";
        }

        if (m == 0) {
            if (t == 3) {
                t = 0;
                s++;

                std::cout << std::endl << "Particle " << s << ": ";
            }
        } else {
            if (t == 4 && s < 4) {
                t = 0;
                s++;

                std::cout << std::endl << "Particle " << s << ": ";
            }
        }
    };

    std::cout << std::endl << "----------------------------------------------------" << std::endl;
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

#ifndef NWA
int calc_me_kinematics_from_int(pStat self, double mH2, double Thb1, double Phb1, double Rhb1, double Thb1b, double Phb1b, double Rhb2, double Thb2)
#else
int calc_me_kinematics_from_int(pStat self, double Thb1, double Phb1, double Rhb1, double Thb1b, double Phb1b, double Rhb2, double Thb2)
#endif
{
    auto p = reinterpret_cast<calc_me*>(self);
    
    #ifndef NWA
        return p->calc_kinematics_from_int(mH2, Thb1, Phb1, Rhb1, Thb1b, Phb1b, Rhb2, Thb2);
    #else
        return p->calc_kinematics_from_int(Thb1, Phb1, Rhb1, Thb1b, Phb1b, Rhb2, Thb2);
    #endif
}

/**
 * matrix_element: [0] -> MG5; [1] -> Physsim
 * 
 * 
*/
double* calc_me_mc_batch(pStat self, double reco_kin[], double int_variables[], int n_elements, int use_transer_funcs, int me_type)
{
    double *buffer;
    buffer = (double*)malloc(n_elements * sizeof(double));

    auto p = reinterpret_cast<calc_me*>(self);
    
    #ifdef ZHH
    bool is_signal = true;
    #else
    bool is_signal = false;
    #endif

    p->me_set(me_type, is_signal);
    p->mc_batch(reco_kin, int_variables, n_elements, buffer, use_transer_funcs);

    return buffer;
}

void calc_me_mem_init(pStat self, double evt_constants[]){
    auto p = reinterpret_cast<calc_me*>(self);
    return p->mem_init(evt_constants);
}

void calc_me_me_set(pStat self, int me_type, bool is_signal){
    auto p = reinterpret_cast<calc_me*>(self);
    return p->me_set(me_type, is_signal);
}

double calc_me_rambo(pStat self)
{
    auto p = reinterpret_cast<calc_me*>(self);
    return p->calc_me_rambo();
}
