//lcme::LCMEZHH *_zhh;

#include "CalculateME.hpp"

#include <iostream>
#include <vector>
#include <string>
#include <cmath>

//----------- Implementation of class lcmezhh -----------//

calc_zhh::calc_zhh(double H_mass, double pol_e, double pol_p, int z_decay_mode, int me_type)
{
    _zhh = new lcme::LCMEZHH("LCMEZHH", "ZHH", H_mass, pol_e, pol_p);
    _zhh->SetZDecayMode(z_decay_mode); // to muon; see LCME code for internal mappings to PDG ids

    _zhh->SetNoHDecay(2);
    _zhh->SetPropagator(1);
    _zhh->SetMEType(me_type);


    std::cout << " [TRACE] LCMEZHH created Ok." << std::endl;
}

calc_zhh::~calc_zhh()
{
    std::cout << " [TRACE] LCMEZHH deleted OK" << std::endl;
}

// Px,Py,Pz,E
double calc_zhh::calc_zhh_calc(double l1_fm[], double l2_fm[], double H1_fm[], double H2_fm[]) const
{
    TLorentzVector zhh_lortz[4] = {
        TLorentzVector(l1_fm),
        TLorentzVector(l2_fm),
        TLorentzVector(H1_fm),
        TLorentzVector(H2_fm)
    };

    _zhh->SetMomentumFinal(zhh_lortz);

    double result = _zhh->GetMatrixElement2();

    std::cout << result << std::endl;

    return result;
}

//----------- Implementation of class calc_zzh -----------//

calc_zzh::calc_zzh(double H_mass, double pol_e, double pol_p, int z1_decay_mode, int z2_decay_mode, int me_type)
{
    _zzh = new lcme::LCMEZZH("LCMEZZH", "ZZH", H_mass, pol_e, pol_p);
    _zzh->SetZDecayMode(z1_decay_mode, z2_decay_mode);

    _zzh->SetNoZDecay(0);
    _zzh->SetNoHDecay(1);
    _zzh->SetPropagator(1);
    _zzh->SetMEType(me_type);

    std::cout << " [TRACE] LCMEZZH created Ok." << std::endl;
}

calc_zzh::~calc_zzh()
{
    std::cout << " [TRACE] LCMEZZH deleted OK" << std::endl;
}

double calc_zzh::calc_zzh_calc(double l1_fm[], double l2_fm[], double z2_d1_fm[], double z2_d2_fm[], double H_fm[]) const
{
    TLorentzVector zzh_lortz[5] = {
        TLorentzVector(l1_fm),
        TLorentzVector(l2_fm),
        TLorentzVector(z2_d1_fm),
        TLorentzVector(z2_d2_fm),
        TLorentzVector(H_fm)
    };

    _zzh->SetMomentumFinal(zzh_lortz);

    double result = _zzh->GetMatrixElement2();

    std::cout << result << std::endl;

    return result;
}


//---------- C-Interface for class lcmezhh ---------------------//

pStat calc_zhh_new(double H_mass, double pol_e, double pol_p, int z_decay_mode, int me_type)
// METype = 2: matrix element only
// METype = 1: differential cross section
{
    return new (std::nothrow) calc_zhh(H_mass, pol_e, pol_p, z_decay_mode, me_type);
}

void calc_zhh_del(pStat self)
{
    delete reinterpret_cast<calc_zhh*>(self);
}

double calc_zhh_calc(pStat self, double l1_fm[], double l2_fm[], double H1_fm[], double H2_fm[])
{
    auto p = reinterpret_cast<calc_zhh*>(self);
    return p->calc_zhh_calc(l1_fm, l2_fm, H1_fm, H2_fm);
}

//---------- C-Interface for class lcmezzh ---------------------//

pStat calc_zzh_new(double H_mass, double pol_e, double pol_p, int z1_decay_mode, int z2_decay_mode, int me_type)
{
    return new (std::nothrow) calc_zzh(H_mass, pol_e, pol_p, z1_decay_mode, z2_decay_mode, me_type);
}

void calc_zzh_del(pStat self)
{
    delete reinterpret_cast<calc_zzh*>(self);
}

double calc_zzh_calc(pStat self, double l1_fm[], double l2_fm[], double z2_d1_fm[], double z2_d2_fm[], double H_fm[])
{
    auto p = reinterpret_cast<calc_zzh*>(self);
    return p->calc_zzh_calc(l1_fm, l2_fm, z2_d1_fm, z2_d2_fm, H_fm);
}


//---------- Helper Functions ---------------------//

// Only valid for e+e- @ 500 GeV; see Physsim
int getZDecayModeFromPDG(int pdg)
{
  int pdgs[11] = {12, 14, 16, 11, 13, 15, 2, 4, 1, 3 ,  5};
  int zdms[11] = {1 , 2,  3,  4,  5,  6,  7, 8, 9, 10, 11};

  for (unsigned int i = 0; i < 11; i++) {
    if (pdgs[i] == pdg) {
      return zdms[i];
    }
  }

  return 0;
}