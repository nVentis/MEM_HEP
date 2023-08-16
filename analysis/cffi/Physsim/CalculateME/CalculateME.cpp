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


    //std::cout << " [TRACE] LCMEZHH created Ok." << std::endl;
}

calc_zhh::~calc_zhh()
{
    //std::cout << " [TRACE] LCMEZHH deleted OK" << std::endl;
}

/**
 * E,Px,Py,Pz
*/
double calc_zhh::calc_zhh_calc(double momenta[]) const
{
    TLorentzVector zhh_lortz[5] = {
        TLorentzVector(momenta[ 1], momenta[ 2], momenta[ 3], momenta[ 0]),
        TLorentzVector(momenta[ 5], momenta[ 6], momenta[ 7], momenta[ 4]),
        TLorentzVector(momenta[ 9], momenta[10], momenta[11], momenta[ 8]),
        TLorentzVector(momenta[13], momenta[14], momenta[15], momenta[12]),
    };

    _zhh->SetMomentumFinal(zhh_lortz);

    return _zhh->GetMatrixElement2();
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

    //std::cout << " [TRACE] LCMEZZH created Ok." << std::endl;
}

calc_zzh::~calc_zzh()
{
    //std::cout << " [TRACE] LCMEZZH deleted OK" << std::endl;
}

/**
 * momenta[]: for lep1, lep2, z2d1, z2d2, H
 * order: E,x,y,z
 * 
 * 
*/
double calc_zzh::calc_zzh_calc(double momenta[]) const
{
    TLorentzVector zzh_lortz[5] = {
        TLorentzVector(momenta[ 1], momenta[ 2], momenta[ 3], momenta[ 0]),
        TLorentzVector(momenta[ 5], momenta[ 6], momenta[ 7], momenta[ 4]),
        TLorentzVector(momenta[ 9], momenta[10], momenta[11], momenta[ 8]),
        TLorentzVector(momenta[13], momenta[14], momenta[15], momenta[12]),
        TLorentzVector(momenta[17], momenta[18], momenta[19], momenta[16])
    };

    _zzh->SetMomentumFinal(zzh_lortz);

    return _zzh->GetMatrixElement2();
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

/**
 * l1_fm[4], l2_fm[4], H1_fm[4], H2_fm[4]
 * -> double momenta[16]
*/
double calc_zhh_calc(pStat self, double momenta[])
{
    auto p = reinterpret_cast<calc_zhh*>(self);
    return p->calc_zhh_calc(momenta);
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

/**
 * l1_fm[4], l2_fm[4], z2_d1_fm[4], z2_d2_fm[4], H_fm[4]
 * -> double momenta[20]
*/
double calc_zzh_calc(pStat self, double momenta[])
{
    auto p = reinterpret_cast<calc_zzh*>(self);
    return p->calc_zzh_calc(momenta);
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