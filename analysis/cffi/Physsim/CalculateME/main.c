#include "CalculateME.hpp"

double calc_zhh_single(double momenta[])
{
    // Internal
    pStat calc = calc_zhh_new(125., -1., 1, 5, 2);

    double me_zhh = calc_zhh_calc(calc, momenta);

    calc_zhh_del(calc);

    return me_zhh;
}

double calc_zzh_single(double momenta[])
{
    pStat calc = calc_zzh_new(125., -1., 1., 5, 11, 2);

    double me = calc_zzh_calc(calc, momenta);

    calc_zzh_del(calc);

    return me;
}