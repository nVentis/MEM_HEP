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

#include "physsim/LCMEZHH.h"
#include "physsim/LCMEZZH.h"
#include "TLorentzVector.h"

class calc_zhh {
  lcme::LCMEZHH *_zhh;

  public:
    calc_zhh(double H_mass, double pol_e, double pol_p, int z_decay_mode, int me_type);
    calc_zhh(calc_zhh const&)            = delete;
    calc_zhh& operator=(calc_zhh const&) = delete;
    ~calc_zhh();

    double calc_zhh_calc(double momenta[]) const;

};

class calc_zzh {
  lcme::LCMEZZH *_zzh;

  public:
    calc_zzh(double H_mass, double pol_e, double pol_p, int z1_decay_mode, int z2_decay_mode, int me_type);
    calc_zzh(calc_zzh const&)            = delete;
    calc_zzh& operator=(calc_zzh const&) = delete;
    ~calc_zzh();

    double calc_zzh_calc(double momenta[]) const;

};

#endif //-- End of __cplusplus definition //

//============ C-interface for class calc_zhh ============//

// Opaque pointer type alias for C-lang
typedef void* pStat;

EXPORT_C pStat   calc_zhh_new(double H_mass, double pol_e, double pol_p, int z_decay_mode, int me_type);
EXPORT_C void    calc_zhh_del (pStat self);
EXPORT_C double  calc_zhh_calc(pStat self, double momenta[]);

//============ C-interface for class calc_zzh ============//

EXPORT_C pStat   calc_zzh_new(double H_mass, double pol_e, double pol_p, int z1_decay_mode, int z2_decay_mode, int me_type);
EXPORT_C void    calc_zzh_del (pStat self);
EXPORT_C double  calc_zzh_calc(pStat self, double momenta[]);


#endif