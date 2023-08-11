#include "CalcME.hpp"

double calc(int iters)
{
    /*
    // Px,Py,Pz,E; as from TLorentzVector
    // Input
    double l1_fm[] = {2., 3., 4., 55.};
    double l2_fm[] = {3., 4., 5., 26.};

    double dec1_fm[] = {4., 5., 6., 37.};
    double dec2_fm[] = {5., 6., 7., 48.};

    double dec3_fm[] = {6., 7., 8., 59.};
    double dec4_fm[] = {7., 8., 9., 48.};
    */

    // Internal
    pStat calc_me = calc_me_new();

    /*
    double Bos_fm[4] = {
        dec1_fm[0] + dec2_fm[0],
        dec1_fm[1] + dec2_fm[1],
        dec1_fm[2] + dec2_fm[2],
        dec1_fm[3] + dec2_fm[3]
    };

    double me_zhh = calc_zhh_calc(calc_zhh, l1_fm, l2_fm, Bos_fm, H_fm);
    */

    double some_me = 0.;

    for (int i = 0; i < iters; i++) {
        some_me += calc_me_rambo(calc_me);
    }
    
    calc_me_del(calc_me);

    return some_me;
}