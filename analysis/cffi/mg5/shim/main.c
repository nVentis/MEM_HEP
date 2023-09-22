#include "CalcME.hpp"

/**
 * param_card: path to a param card, use /afs/desy.de/user/b/bliewert/public/MarlinWorkdirs/MEM_HEP/analysis/cffi/mg5/mg5/Cards/param_card.dat
 * helicity_selection int[2n]: with particle @(2n) and helicity @(2n+1)
 * selected_helicities int: len(helicity_selection)/2
*/
double calc_rambo(const char param_card[], double energy, int helicity_selection[], int selected_helicities)
{
    // Internal
    pStat calc_me = calc_me_new(param_card, energy);

    if (sizeof(helicity_selection) > 0) {
        for (int j = 0; j < selected_helicities; j++) {
            calc_me_set_helicity(calc_me, helicity_selection[2*j], helicity_selection[2*j+1]);
        }
    }

    double some_me = calc_me_rambo(calc_me);
    
    calc_me_del(calc_me);

    return some_me;
}

/**
 * param_card: path to a param card, use /afs/desy.de/user/b/bliewert/public/MarlinWorkdirs/MEM_HEP/analysis/cffi/mg5/mg5/Cards/param_card.dat
 * helicity_selection int[2n]: with particle @(2n) and helicity @(2n+1)
 * selected_helicities int: len(helicity_selection)/2
 * momenta double[24*l]: flat final state momenta for l events; for each event order is mu-, mu+, b1, bbar1, b2, bbar2 with each E,Px,Py,Pz
 * n_elements: length of result, i.e. double[n_elements]; should be len(momenta)/24
*/
double* calc(const char param_card[], double energy, int helicity_selection[], int selected_helicities, double momenta[], int n_elements)
{
    // Internal
    pStat calc_me = calc_me_new(param_card, energy);

    if (sizeof(helicity_selection) > 0) {
        for (int j = 0; j < selected_helicities; j++) {
            calc_me_set_helicity(calc_me, helicity_selection[2*j], helicity_selection[2*j+1]);
        }
    }

    double *some_me = calc_me_multi(calc_me, momenta, n_elements);
    
    calc_me_del(calc_me);

    return some_me;
}

/**
 * param_card: path to a param card, use /afs/desy.de/user/b/bliewert/public/MarlinWorkdirs/MEM_HEP/analysis/cffi/mg5/mg5/Cards/param_card.dat
 * helicity_selection int[2n]: with particle @(2n) and helicity @(2n+1)
 * selected_helicities int: len(helicity_selection)/2
 * momenta double[24*l]: flat final state momenta for l events; for each event order is mu-, mu+, b1, bbar1, b2, bbar2 with each E,Px,Py,Pz
 * n_elements: length of result, i.e. double[n_elements]; should be len(momenta)/24
*/
double* calc_mc_batch(const char param_card[], double energy, int helicity_selection[], int selected_helicities, double reco_kin[], double int_variables[], int n_elements)
{
    // Internal
    pStat calc_me = calc_me_new(param_card, energy);

    if (sizeof(helicity_selection) > 0) {
        for (int j = 0; j < selected_helicities; j++) {
            calc_me_set_helicity(calc_me, helicity_selection[2*j], helicity_selection[2*j+1]);
        }
    }

    double *some_me = calc_me_mc_batch(calc_me, reco_kin, int_variables, n_elements);
    
    calc_me_del(calc_me);

    return some_me;
}