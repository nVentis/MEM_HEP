#ifndef _CODECPP_
#define _CODECPP_

#ifdef __cplusplus
  #define EXPORT_C extern "C"
#else
  #define EXPORT_C
  #include <stddef.h>
#endif

//============ C++ Only Header =================//
#ifdef __cplusplus  // Enabled only for C++ compilers
#define _USE_MATH_DEFINES
#include <iostream>

#include "CPPProcess.h"
#include "rambo.h"

#include <ios>
#include <iomanip>
#include <vector>
#include <array>
#include <string>
#include <cmath>
#include <memory>

// MEElement
#include "physsim/LCMEZHH.h"
#include "physsim/LCMEZZH.h"
#include "TLorentzVector.h"

struct PhysConstants {
  double mb = 4.8;
  double epsilon = 1e-3;
  double dEmax = .1;
  double dpmax = .1;
  double sqrt_s = 500.;
  double system_E = 0.;
  double system_px = 0.;
  double system_py = 0.;
  double system_pz = 0.;
  double mH = 125.;
};

/*
class MEElement
{
    protected:
        PhysConstants *constants;
        
    public:
        virtual ~MEElement() {};
        virtual double single(std::vector<double*> momenta);

        void setConstants(PhysConstants *consts) { constants = consts; };
};
*/

class MEElementPhyssimSig
{
    protected:
        PhysConstants *constants;

    private:
        lcme::LCMEZHH *_lcme = nullptr;

    public:
        MEElementPhyssimSig(PhysConstants *consts, double pol_e, double pol_p, int z_decay_mode);
        ~MEElementPhyssimSig() { delete _lcme; };

        double single(std::vector<double*> momenta);
        void setConstants(PhysConstants *consts) { constants = consts; };
};

class MEElementPhyssimBkg
{
    protected:
        PhysConstants *constants;

    private:
        lcme::LCMEZZH *_lcme = nullptr;

    public:
        MEElementPhyssimBkg(PhysConstants *consts, double pol_e, double pol_p, int z1_decay_mode, int z2_decay_mode);
        ~MEElementPhyssimBkg() { delete _lcme; };

        double single(std::vector<double*> momenta);
        void setConstants(PhysConstants *consts) { constants = consts; };
};












using vec3 = std::array<double,3>;
using kinematics = std::array<double, 28>; // 4*3 spherical coordinates, then 4*4 four vectors

std::array<double, 2> get_angles(double x, double y, double z){
  if (x == 0. && y == 0. && z == 0.) {
    // to handle use_transer_funcs=false case
    return std::array<double, 2> {
      0,
      0
    };
  } else {
    return std::array<double, 2> {
      std::acos(z/std::sqrt(std::pow(x, 2.) + std::pow(y, 2.) + std::pow(z, 2.))),
      std::atan2(y, x)
    };
  }
}

std::array<std::array<double, 2>, 4> get_reco_angles(double reco_kin[]) {
  std::array<std::array<double, 2>, 4> res;
  // reco_kin[0-3]: mu-
  // reco_kin[4-7]: mu+

  int i;
  for (i = 0; i < 4; i++)
    res[i] = get_angles(reco_kin[9+(i*4)], reco_kin[10+(i*4)], reco_kin[11+(i*4)]);

  return res;
}

const vec3 vec3_sph(double theta, double phi, double len) {
  return vec3 {
    len*std::sin(theta)*std::cos(phi),
    len*std::sin(theta)*std::sin(phi),
    len*std::cos(theta)
  };
};

vec3 vec3_add(vec3 a, vec3 b) {
  return vec3{
    a[0] + b[0],
    a[1] + b[1],
    a[2] + b[2]
  };
};

const double vec3_dot(vec3 a, vec3 b) {
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

class calc_me {
  private:
    CPPProcess *_cppp;

    // Transfer function related
    // default values estimated from MC samples
    // BKGHYP and SIGHYP to allow separate transfer functions for both hypotheses (TESTING)
    #ifdef BKGHYP
    double tf_E1_args[2] { -1.27071409, 6.54144637 }; // H
    double tf_E2_args[2] { -0.14618351, 5.41285247 }; // Z

    double tf_Th1_args[2] { -5.76347284e-05,  2.02238842e-02 };
    double tf_Th2_args[2] { 6.79667849e-06, 1.91973789e-02 };
    
    double tf_Ph1_args[2] { -0.00030157,  0.02824722 };
    double tf_Ph2_args[2] { 0.00090488, 0.02750561 };

    //double tf_E_args[2] { -0.60700026, 5.89539517};
    //double tf_Th_args[2] { -6.19007451e-06, 1.97328833e-02 };
    //double tf_Ph_args[2] { 5.19239786e-05, 2.83442591e-02 };
    #else
      #ifdef SIGHYP
        double tf_E1_args[2] { -1.31771938, 6.44408755 }; // H
        double tf_E2_args[2] { -1.35264493, 6.41472459 }; // H

        double tf_Th1_args[2] { 3.14744981e-06, 1.95730114e-02 };
        double tf_Th2_args[2] { -7.81672016e-05,  1.95288882e-02 };

        double tf_Ph1_args[2] { 0.00024616, 0.02817029 };
        double tf_Ph2_args[2] { 0.00082556, 0.02781802 };

        //double tf_E_args[2] { -1.29695176,  6.34412724 };
        //double tf_Th_args[2] { -4.03399692e-05, 1.96212741e-02 };
        //double tf_Ph_args[2] { 0.00039051, 0.02811957 };
      #else
        double tf_E1_args[2] { -1.12594285, 6.24937028};
        double tf_E2_args[2] { -1.12594285, 6.24937028};

        double tf_Th1_args[2] { -3.90908961e-05, 1.96662831e-02 };
        double tf_Th2_args[2] { -3.90908961e-05, 1.96662831e-02 };

        double tf_Ph1_args[2] { 0.0001748, 0.02819419 };
        double tf_Ph2_args[2] { 0.0001748, 0.02819419 };

        //double tf_E_args[2] { -1.12594285, 6.24937028};
        //double tf_Th_args[2] { -3.90908961e-05, 1.96662831e-02 };
        //double tf_Ph_args[2] { 0.0001748, 0.02819419 };
      #endif
    #endif

    PhysConstants constants; // TODO

    // MEM related
    double calc_jac(double Rhb1, double Thb1, double Phb1,
                    double Rhb1b, double Thb1b, double Phb1b,
                    double Rhb2, double Thb2, double Phb2,
                    double Rhb2b, double Thb2b, double Phb2b);
    kinematics kin{};

    int err_map[11]{ -1, -2, -3, -4, -5, -6, -7, -8, -9, -10 , -11 };

    double calc_tf_E(int jet_idx, double a, double b);
    double calc_tf_Th(int jet_idx, double a, double b);
    double calc_tf_Ph(int jet_idx, double a, double b);

    void kin_debug_print();

    // Other often needed stuff
    double mb_pow2 = std::pow(4.8, 2.);
    bool mem_init_run{false};

  protected:
    #ifdef ZHH
    MEElementPhyssimSig *me_element = nullptr;
    #else
    MEElementPhyssimBkg *me_element = nullptr;
    #endif

  public:
    calc_me(std::string param_card, double energy);
    calc_me(calc_me const&)            = delete;
    calc_me& operator=(calc_me const&) = delete;
    ~calc_me() {
      //delete me_element; // some errors; TODO
      delete _cppp;
    };

    void set_helicity(int particle, int helicity);
    double calc_me_single_full_kin(std::vector<double*> momenta);
    double* calc_me_multi(double momenta[], int n_elements, double buffer[]);
    double calc_me_rambo();

    // MEM related
    void mem_init(double evt_constants[]);
    void me_set(int me_type, bool is_signal); // 0: MG5; 1: Physsim
    kinematics get_kinematics() { return kin; };

    void mc_batch(double reco_kin[], double int_variables[], int n_elements, double buffer[], int use_transer_funcs);

    // Public mainly for testing; kin contains main information
    #ifndef NWA
    int calc_kinematics_from_int(double mH2, double Thb1, double Phb1, double Rhb1, double Thb1b, double Phb1b, double Rhb2, double Thb2);
    #else
    int calc_kinematics_from_int(double Thb1, double Phb1, double Rhb1, double Thb1b, double Phb1b, double Rhb2, double Thb2);
    #endif
};

#endif //-- End of __cplusplus definition //

//============ C-interface for class calc_me ============//

// Opaque pointer type alias for C-lang
typedef void* pStat;

EXPORT_C pStat   calc_me_new                  (const char param_card[], double energy);
EXPORT_C void    calc_me_set_helicity         (pStat self, int particle, int helicity);
EXPORT_C void    calc_me_del                  (pStat self);
EXPORT_C double* calc_me_multi                (pStat self, double momenta[], int n_elements);
#ifndef NWA
EXPORT_C int     calc_me_kinematics_from_int  (pStat self, double mH2, double Thb1, double Phb1, double Rhb1, double Thb1b, double Phb1b, double Rhb2, double Thb2);
#else
EXPORT_C int     calc_me_kinematics_from_int  (pStat self, double Thb1, double Phb1, double Rhb1, double Thb1b, double Phb1b, double Rhb2, double Thb2);
#endif
EXPORT_C double* calc_me_mc_batch             (pStat self, double reco_kin[], double int_variables[], int n_elements, int use_transer_funcs, int me_type);
EXPORT_C void    calc_me_mem_init             (pStat self, double evt_constants[]);
EXPORT_C void    calc_me_me_set               (pStat self, int me_type);
EXPORT_C double  calc_me_rambo                (pStat self);

#endif