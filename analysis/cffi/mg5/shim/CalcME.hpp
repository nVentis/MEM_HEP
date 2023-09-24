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

struct phys_const {
  double mb = 4.8;
  double epsilon = 1e-3;
  double dEmax = .1;
  double dpmax = .1;
  double sqrt_s = 500.;
  double system_E = 0.;
  double system_px = 0.;
  double system_py = 0.;
  double system_pz = 0.;
};

using vec3 = std::array<double,3>;
using kinematics = std::array<double, 28>; // 4*3 spherical coordinates, then 4*4 four vectors

std::array<double, 2> get_angles(double x, double y, double z){
  std::array<double, 2> res {
    std::acos(z/std::sqrt(std::pow(x, 2.) + std::pow(y, 2.) + std::pow(z, 2.))),
    std::atan2(y, x)
  };
  
  return res;
}

std::array<std::array<double, 2>, 4> get_reco_angles(double reco_kin[]) {
  std::array<std::array<double, 2>, 4> res;
  // reco_kin[0-3]: mu-
  // reco_kin[4-7]: mu+

  int i;
  for (i = 0; i < 4; i++)
    res[i] = get_angles(reco_kin[8+(i*4)], reco_kin[9+(i*4)], reco_kin[10+(i*4)]);

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
    double tf_E_args[2] { -1.12594285, 6.24937028};
    double tf_Th_args[2] { -3.90908961e-05, 1.96662831e-02 };
    double tf_Ph_args[2] { 0.0001748, 0.02819419 };

    phys_const constants{}; // TODO

    // MEM related
    double calc_jac(double Rhb1, double Thb1, double Phb1,
                    double Rhb1b, double Thb1b, double Phb1b,
                    double Rhb2, double Thb2, double Phb2,
                    double Rhb2b, double Thb2b, double Phb2b);
    kinematics kin{};

    int err_map[11]{ -1, -2, -3, -4, -5, -6, -7, -8, -9, -10 , -11 };

    double calc_tf_E(double a, double b);
    double calc_tf_Th(double a, double b);
    double calc_tf_Ph(double a, double b);

    void kin_debug_print();

    // Other often needed stuff
    double mb_pow2 = std::pow(4.8, 2.);

  public:
    calc_me(std::string param_card, double energy);
    calc_me(calc_me const&)            = delete;
    calc_me& operator=(calc_me const&) = delete;
    ~calc_me();

    void set_helicity(int particle, int helicity);
    double calc_me_single_full_kin(std::vector<double*> momenta);
    double* calc_me_multi(double momenta[], int n_elements, double buffer[]);
    double calc_me_rambo();

    // MEM related
    void mem_init(double evt_constants[]);
    kinematics get_kinematics() { return kin; };

    void mc_batch(double reco_kin[], double int_variables[], int n_elements, double buffer[]);

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
EXPORT_C double* calc_me_mc_batch             (pStat self, double reco_kin[], double int_variables[], int n_elements);
EXPORT_C void    calc_me_mem_init             (pStat self, double evt_constants[]);
EXPORT_C double  calc_me_rambo                (pStat self);

#endif