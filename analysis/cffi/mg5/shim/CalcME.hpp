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
  double mb;
  double sqrt_s;
  double system_px:
  double system_py;
  double system_pz;
};

using vec3 = std::array<double,3>;

const std::array<double, 2> get_angles(double x, double y, double z) {
}

const vec3 spherical_vec(double theta, double phi, double len) {
  return vec3 {
    len*std::sin(theta)*std::cos(phi),
    len*std::sin(theta)*std::sin(phi),
    len*std::cos(phi)
  };
};

const double vec3_dot(vec3 a, vec3 b) {
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

class calc_me {
  private:
    CPPProcess *_cppp;
    double E{};

    // Transfer function related
    double tf_E_args[2];
    double tf_Th_args[2];
    double tf_Ph_args[2];

    phys_const constants { 4.8 }; // TOTO: setConstants

    // MEM related
    double calc_jac(double params[]);
    double get_kinematics_from_int(double int_variables[]);

    int err_map[11]{ -1, -2, -3, -4, -5, -6, -7, -8, -9, -10 , -11 };

    // Integration variables
    // TODO: narrow width approximation to get rid of mH2 integration for ZHH and ZZH; however need to adjust MadGraph matrix element
    double mH2;
    double Thb1;
    double Phb1;
    double Rhb1;
    double Thb1b;
    double Phb1b;
    double Rhb2;
    double Thb2;

    // Other often needed stuff
    double mb_pow2 = std::pow(4.8, 2.);

  public:
    calc_me(std::string param_card, double energy);
    calc_me(calc_me const&)            = delete;
    calc_me& operator=(calc_me const&) = delete;
    ~calc_me();

    void set_helicity(int particle, int helicity);
    double* calc_me_multi(double momenta[], int n_elements, double buffer[]);
    double calc_me_rambo();

};

#endif //-- End of __cplusplus definition //

//============ C-interface for class calc_me ============//

// Opaque pointer type alias for C-lang
typedef void* pStat;

EXPORT_C pStat   calc_me_new          (const char param_card[], double energy);
EXPORT_C void    calc_me_set_helicity (pStat self, int particle, int helicity);
EXPORT_C void    calc_me_del          (pStat self);
EXPORT_C double* calc_me_multi        (pStat self, double momenta[], int n_elements);
EXPORT_C double  calc_me_rambo        (pStat self);

#endif