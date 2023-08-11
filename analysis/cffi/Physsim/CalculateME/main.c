#include "CalculateME.hpp"

double calc()
{
    // Px,Py,Pz,E; as from TLorentzVector
    // Input
    double l1_fm[] = {2., 3., 4., 55.};
    double l2_fm[] = {3., 4., 5., 26.};

    double dec1_fm[] = {4., 5., 6., 37.};
    double dec2_fm[] = {5., 6., 7., 48.};

    double dec3_fm[] = {6., 7., 8., 59.};
    double dec4_fm[] = {7., 8., 9., 48.};

    // Internal
    pStat calc_zhh = calc_zhh_new(125., 1., 0, 5, 2);
    pStat calc_zzh = calc_zzh_new(125., 1., 0, 5, 4, 2);

    double Bos_fm[4] = {
        dec1_fm[0] + dec2_fm[0],
        dec1_fm[1] + dec2_fm[1],
        dec1_fm[2] + dec2_fm[2],
        dec1_fm[3] + dec2_fm[3]
    };

    double H_fm[4] = {
        dec3_fm[0] + dec4_fm[0],
        dec3_fm[1] + dec4_fm[1],
        dec3_fm[2] + dec4_fm[2],
        dec3_fm[3] + dec4_fm[3]
    };

    double me_zhh = calc_zhh_calc(calc_zhh, l1_fm, l2_fm, Bos_fm, H_fm);
    double me_zzh = calc_zzh_calc(calc_zzh, l1_fm, l2_fm, dec1_fm, dec2_fm, H_fm);

    calc_zhh_del(calc_zhh);
    calc_zzh_del(calc_zzh);

    return me_zhh;
}