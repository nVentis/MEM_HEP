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
#include <string>
#include <cmath>

class calc_me {
  private:
    CPPProcess *_cppp;

  public:
    calc_me();
    calc_me(calc_me const&)            = delete;
    calc_me& operator=(calc_me const&) = delete;
    ~calc_me();

    double calc_me_calc();
    double calc_me_rambo();

};

#endif //-- End of __cplusplus definition //

//============ C-interface for class calc_me ============//

// Opaque pointer type alias for C-lang
typedef void* pStat;

EXPORT_C pStat   calc_me_new  ();
EXPORT_C void    calc_me_del  (pStat self);
EXPORT_C double  calc_me_calc (pStat self);
EXPORT_C double  calc_me_rambo(pStat self);

#endif