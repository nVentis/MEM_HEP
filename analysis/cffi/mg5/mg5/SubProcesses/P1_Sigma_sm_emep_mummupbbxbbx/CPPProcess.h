//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph5_aMC@NLO v. 3.5.1, 2023-07-11
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#ifndef MG5_Sigma_sm_emep_mummupbbxbbx_H
#define MG5_Sigma_sm_emep_mummupbbxbbx_H

#include <complex> 
#include <vector> 

#include "Parameters_sm.h"

using namespace std; 

//==========================================================================
// A class for calculating the matrix elements for
// Process: e- e+ > mu- mu+ h h WEIGHTED<=8 @1
// *   Decay: h > b b~ WEIGHTED<=2
// *   Decay: h > b b~ WEIGHTED<=2
//--------------------------------------------------------------------------

class CPPProcess
{
  public:

    // Constructor.
    CPPProcess() {}

    // Initialize process.
    virtual void initProc(string param_card_name); 

    // Calculate flavour-independent parts of cross section.
    virtual void sigmaKin(); 

    // Evaluate sigmaHat(sHat).
    virtual double sigmaHat(); 

    // Info on the subprocess.
    virtual string name() const {return "e- e+ > mu- mu+ b b~ b b~ (sm)";}

    virtual int code() const {return 1;}

    const vector<double> & getMasses() const {return mME;}

    // Get and set momenta for matrix element evaluation
    vector < double * > getMomenta(){return p;}
    void setMomenta(vector < double * > & momenta){p = momenta;}
    void setInitial(int inid1, int inid2){id1 = inid1; id2 = inid2;}

    // Get matrix element vector
    const double * getMatrixElements() const {return matrix_element;}

    // Constants for array limits
    static const int ninitial = 2; 
    static const int nexternal = 8; 
    static const int nprocesses = 1; 

    // MANUAL HELICITY SELECTION START
    void setHelicity(int particle, int helicity){ selected_helicities[particle] = helicity; }
    int checkHelicityComb(const int helicities[]) {
      for (int j = 0; j < nexternal; j++) {
        if (selected_helicities[j] != 0 && selected_helicities[j] != helicities[j]) return 0;
      };
      return 1;
    };
    int selected_helicities[nexternal] = {nexternal * 0};
    // MANUAL HELICITY SELECTION END

    // MANUAL LAMBDA START
    void setLambda(double lambda) {
      // default: mdl_lam = mdl_MH__exp__2/(2. * mdl_vev__exp__2); 
      pars->mdl_lam = lambda;
      pars->mdl_muH = sqrt(pars->mdl_lam * pars->mdl_vev__exp__2);
      pars->GC_69 = -6. * pars->mdl_complexi * pars->mdl_lam * pars->mdl_vev; // Higgs self-coupling
    };
    // MANUAL LAMBDA END

  private:

    // Private functions to calculate the matrix element for all subprocesses
    // Calculate wavefunctions
    void calculate_wavefunctions(const int perm[], const int hel[]); 
    static const int nwavefuncs = 15; 
    std::complex<double> w[nwavefuncs][18]; 
    static const int namplitudes = 4; 
    std::complex<double> amp[namplitudes]; 
    double matrix_1_emep_mummuphh_h_bbx_h_bbx(); 

    // Store the matrix element value from sigmaKin
    double matrix_element[nprocesses]; 

    // Color flows, used when selecting color
    double * jamp2[nprocesses]; 

    // Pointer to the model parameters
    Parameters_sm * pars; 

    // vector with external particle masses
    vector<double> mME; 

    // vector with momenta (to be changed each event)
    vector < double * > p; 
    // Initial particle ids
    int id1, id2;     

}; 


#endif  // MG5_Sigma_sm_emep_mummupbbxbbx_H
