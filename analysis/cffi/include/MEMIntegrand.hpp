#include "vector"
#include "array"
#include "cmath"
#include "iostream"
#include "PhysConstants.hpp"
#include "MEElement.hpp"
#include "vec3.hpp"
#include "MEMTransferFunction.hpp"

using kinematics = std::array<double, 28>; // 4*3 spherical coordinates, then 4*4 four vectors

class MEMIntegrand
{
    private:
        /* data */
        PhysConstants constants;
        MEElement melement;
        double mb_pow2;
        kinematics kin{};

        // Transfer function parameterization (hard-coded: Lorentzian); default values estimated from MC samples
        // BKGHYP and SIGHYP to allow separate transfer functions for both hypotheses (TESTING)
        #ifdef BKGHYP
        double tf_E_args[2] { -0.60700026, 5.89539517};
        double tf_Th_args[2] { -6.19007451e-06, 1.97328833e-02 };
        double tf_Ph_args[2] { 5.19239786e-05, 2.83442591e-02 };
        #else
        #ifdef SIGHYP
            double tf_E_args[2] { -1.29695176,  6.34412724 };
            double tf_Th_args[2] { -4.03399692e-05, 1.96212741e-02 };
            double tf_Ph_args[2] { 0.00039051, 0.02811957 };
        #else
            double tf_E_args[2] { -1.12594285, 6.24937028};
            double tf_Th_args[2] { -3.90908961e-05, 1.96662831e-02 };
            double tf_Ph_args[2] { 0.0001748, 0.02819419 };
        #endif
        #endif

    protected:
        // Transfer Functions
        //MEMTransferFunction tf_E;
        //MEMTransferFunction tf_Th;
        //MEMTransferFunction tf_Ph;

    public:
        MEMIntegrand(PhysConstants cons, MEElement me);
        ~MEMIntegrand();

        double calc_jac(double Rhb1, double Thb1, double Phb1,
                    double Rhb1b, double Thb1b, double Phb1b,
                    double Rhb2, double Thb2, double Phb2,
                    double Rhb2b, double Thb2b, double Phb2b);

        #ifndef NWA
        int calc_kinematics_from_int(double mH2, double Thb1, double Phb1, double Rhb1, double Thb1b, double Phb1b, double Rhb2, double Thb2);
        #else
        int calc_kinematics_from_int(double Thb1, double Phb1, double Rhb1, double Thb1b, double Phb1b, double Rhb2, double Thb2);
        #endif

        void kin_debug_print();
        void mc_batch(double reco_kin[], double int_variables[], int n_elements, double buffer[], int use_transer_funcs);

        // Future TODO: setTF(MEMTransferFunction)
        double calc_tf_E(double a, double b);
        double calc_tf_Th(double a, double b);
        double calc_tf_Ph(double a, double b);
        
};

MEMIntegrand::MEMIntegrand(
    PhysConstants cons,
    MEElement me
){
    constants = cons;
    melement = me;
    mb_pow2 = std::pow(constants.mb, 2.);
}

MEMIntegrand::~MEMIntegrand()
{
}
