/**
 * @author Bryan Bliewert (TUM/DESY) <bryan.bliewert@tum.de>
 * @date 20.05.2023
 * 
 * Calculate Matrix Elements (MEs) using ZHH and ZZH processes based on supplied 4-vectors
 * 
 * Parameters:
 * - LepPairCollection<RECONSTRUCTEDPARTICLE>
 * 
 * - HiggsPairCollection<RECONSTRUCTEDPARTICLE>
 * - PreSelectionCollection<RECONSTRUCTEDPARTICLE> [why is this a particle and not a simple boolean flag?]
 * - MCTrueCollection<MCPARTICLE>
 * - Z1DecayMode<int(5)>; defaults to µ+µ- (see ZHH decay modes when running the MEM processor)
 * - HiggsMass<float(125.)>
 * - outputFilename<string("output.root")>
 * - outputTree<string("dataTree")>
*/

#include "CompareMEProcessor.h"
#include <iostream>
#include <vector>
#include <string>
#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>
#include <IMPL/ReconstructedParticleImpl.h>
#include <UTIL/PIDHandler.h>

using namespace lcio ;
using namespace marlin ;
using namespace std ;
using namespace lcme;

template<class T>
TLorentzVector v4(T* p){
  return TLorentzVector( p->getMomentum(), p->getEnergy() );
}

template<class T>
double inv_mass(T* p1, T* p2){
  double e = p1->getEnergy()+p2->getEnergy() ;
  double px = p1->getMomentum()[0]+p2->getMomentum()[0];
  double py = p1->getMomentum()[1]+p2->getMomentum()[1];
  double pz = p1->getMomentum()[2]+p2->getMomentum()[2];
  return( sqrt( e*e - px*px - py*py - pz*pz  ) );
}

// Helper class to delete TrueJet_Parser instance when it comes out of scope
struct DelMe {
  DelMe( std::function<void()> func ) : _func(func) {}
  ~DelMe() { _func(); }
  std::function<void()>  _func;
};

CompareMEProcessor aCompareMEProcessor ;

CompareMEProcessor::CompareMEProcessor() :

  Processor("CompareMEProcessor"),
  m_nRun(0),
  m_nEvt(0),
  m_mode_me(1), // 0-> dsigma, 1-> ME^2
  m_zzh_no_z_decay(0) // 0-> both Zs decay, 1-> 1 Z does not decay
{

	_description = "CompareMEProcessor writes relevant observables to root-file " ;

  registerInputCollection(LCIO::MCPARTICLE,
				 "InputMCTrueCollection",
				 "preselection collection",
				 m_inputMCTrueCollection,
				 std::string("MCParticlesSkimmed")
				 );

  registerInputCollection(LCIO::MCPARTICLE,
				 "InputPreSelection",
				 "preselection collection",
				 m_inputPreSelectionCollection,
				 std::string("preselection")
				 );

  registerInputCollection(LCIO::MCPARTICLE,
				 "InputHdecayMode",
				 "HdecayMode collection with important parameters (flavour)",
				 m_inputHdecayModeCollection,
				 std::string("HdecayMode")
				 );

  registerInputCollection(LCIO::RECONSTRUCTEDPARTICLE,
				 "InputJetCollection",
				 "collection of four jets originating from two higgses",
				 m_inputJetCollection,
				 std::string("RefinedJets")
				 );

  registerInputCollection(LCIO::RECONSTRUCTEDPARTICLE,
				 "InputHiggsPairCollection",
				 "higgs pair collection (two particles)",
				 m_inputHiggsPairCollection,
				 std::string("HiggsPair")
				 );
  
  registerInputCollection(LCIO::RECONSTRUCTEDPARTICLE,
				 "InputLepPairCollection",
				 "lepton pair collection (two leptons)",
				 m_inputLepPairCollection,
				 std::string("LeptonPair")
				 );

  registerInputCollection(LCIO::RECONSTRUCTEDPARTICLE,
				 "JetMatchingCollection",
				 "Name of collection holding parameters about reco/true jet matching",
				 m_inputJetMatchingCollection,
				 std::string("JetMatching")
				 );

  registerProcessorParameter("RequirePreselectionPass",
        "whether (1) or not to skip (1) events that did not pass preselection",
        m_require_presel_pass,
        int(0)
        );

	registerProcessorParameter("Z1DecayPDG",
        "MEM processor mode of decay of Z1 (=Z for ZHH, Z1 for ZZH)",
        m_z1_decay_pdg,
        int(13)
        );

  registerProcessorParameter("Mode",
        "mode of usage; 0:MCParticles, 1:InputHiggsPairCollection with InputJetCollection, 1:TrueJet",
        m_mode,
        int(0) // 0: True/MCParticleSkimmed; 1: HiggsPair (e.g. RefinedJets with given jet pairing) 
        );

  registerProcessorParameter("LeptonMode",
        "for Mode=1,2; specifies which lepton data may be used; 0: same as Mode, 1: LeptonPair, 2: TrueJet",
        m_lepton_mode,
        int(0) // 0: True/MCParticleSkimmed; 1: HiggsPair (e.g. RefinedJets with given jet pairing) 
        );

  registerProcessorParameter("SaveInputKinematics",
        "0: no complete kinematics, 1: save x,y,z,E of all inputs",
        m_saveInputKinematics,
        int(0) // 0: True/MCParticleSkimmed; 1: HiggsPair (e.g. RefinedJets with given jet pairing) 
        );

	registerProcessorParameter("HiggsMass",
        "assumed Hmass",
        m_Hmass,
        float(125.)
        );
	
  registerProcessorParameter("outputFilename",
        "name of output root file",
        m_outputFile,
        std::string("compare_out.root")
        );

  registerProcessorParameter("outputTree",
        "name of resulting TTree in the output file",
        m_outputTree,
        std::string("dataTree")
        );

  // TrueJet_Parser parameters
  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
                           "TrueJets" ,
                           "Name of the TrueJetCollection input collection"  ,
                           _trueJetCollectionName ,
                           std::string("TrueJets") ) ;

  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
                           "FinalColourNeutrals" ,
                           "Name of the FinalColourNeutralCollection input collection"  ,
                           _finalColourNeutralCollectionName ,
                           std::string("FinalColourNeutrals") ) ;

  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
                           "InitialColourNeutrals" ,
                           "Name of the InitialColourNeutralCollection input collection"  ,
                           _initialColourNeutralCollectionName ,
                           std::string("InitialColourNeutrals") ) ;

  registerInputCollection( LCIO::LCRELATION,
                            "TrueJetPFOLink" ,
                            "Name of the TrueJetPFOLink input collection"  ,
                            _trueJetPFOLink,
                            std::string("TrueJetPFOLink") ) ;

  registerInputCollection( LCIO::LCRELATION,
                            "TrueJetMCParticleLink" ,
                            "Name of the TrueJetMCParticleLink input collection"  ,
                            _trueJetMCParticleLink,
                            std::string("TrueJetMCParticleLink") ) ;

  registerInputCollection( LCIO::LCRELATION,
                            "FinalElementonLink" ,
                            "Name of the  FinalElementonLink input collection"  ,
                            _finalElementonLink,
                            std::string("FinalElementonLink") ) ;

  registerInputCollection( LCIO::LCRELATION,
                            "InitialElementonLink" ,
                            "Name of the  InitialElementonLink input collection"  ,
                            _initialElementonLink,
                            std::string("InitialElementonLink") ) ;

  registerInputCollection( LCIO::LCRELATION,
                            "FinalColourNeutralLink" ,
                            "Name of the  FinalColourNeutralLink input collection"  ,
                            _finalColourNeutralLink,
                            std::string("FinalColourNeutralLink") ) ;

  registerInputCollection( LCIO::LCRELATION,
                            "InitialColourNeutralLink" ,
                            "Name of the  InitialColourNeutralLink input collection"  ,
                            _initialColourNeutralLink,
                            std::string("InitialColourNeutralLink") ) ;
}

void CompareMEProcessor::init()
{
  streamlog_out(DEBUG) << "   init called  " << std::endl;
  this->Clear();

  // Configuration
  Double_t pol_e = 0.; // this->parameters()->getFloatVal("beamPol1");
  Double_t pol_p = 0; // this->parameters()->getFloatVal("beamPol2");

  // Runtime variables
  m_nRun = 0;
  m_nEvt = 0;

  m_pTFile = new TFile(m_outputFile.c_str(), "recreate");
  m_pTTree = new TTree(m_outputTree.c_str(), m_outputTree.c_str());
  m_pTTree->SetDirectory(m_pTFile);

  m_pTTree->Branch("run", &m_nRun, "run/I");
  m_pTTree->Branch("event", &m_nEvt, "event/I");
  m_pTTree->Branch("error_code", &m_error_code, "error_code/I");
  /**
   * Error code 
   * 
   * 
   * 
  */

  // 1. Meta
  m_pTTree->Branch("is_zhh", &m_is_zhh, "is_zhh/I");
  m_pTTree->Branch("is_zzh", &m_is_zzh, "is_zzh/I");
  m_pTTree->Branch("passed_preselection", &m_passed_preselection, "passed_preselection/I");
  m_pTTree->Branch("true_h1_decay_pdg", &m_true_h1_decay_pdg, "true_h1_decay_pdg/I");
  m_pTTree->Branch("true_h2_decay_pdg", &m_true_h2_decay_pdg, "true_h2_decay_pdg/I");
  m_pTTree->Branch("true_z2_decay_pdg", &m_true_z2_decay_pdg, "true_z2_decay_pdg/I");

  // 2. Event data
  m_pTTree->Branch("h1_decay_pdg", &m_h1_decay_pdg, "h1_decay_pdg/I");
  m_pTTree->Branch("h2_decay_pdg", &m_h2_decay_pdg, "h2_decay_pdg/I");
  m_pTTree->Branch("z2_decay_pdg", &m_z2_decay_pdg, "z2_decay_pdg/I");

  // 2.a ZHH output
  m_pTTree->Branch("zhh_sigma"  , &m_zhh_sigma  , "zhh_sigma/F");
  m_pTTree->Branch("zhh_sigmall", &m_zhh_sigmall, "zhh_sigmall/F");
  m_pTTree->Branch("zhh_sigmalr", &m_zhh_sigmalr, "zhh_sigmalr/F");
  m_pTTree->Branch("zhh_sigmarl", &m_zhh_sigmarl, "zhh_sigmarl/F");
  m_pTTree->Branch("zhh_sigmarr", &m_zhh_sigmarr, "zhh_sigmarr/F");

  m_pTTree->Branch("zhh_mz"  , &m_zhh_mz  , "zhh_mz/F");
  m_pTTree->Branch("zhh_mhh" , &m_zhh_mhh , "zhh_mhh/F");
  m_pTTree->Branch("zhh_mzhh", &m_zhh_mzhh, "zhh_mzhh/F");

  m_pTTree->Branch("zhh_phi" , &m_zhh_phi , "zhh_phi/F");
  m_pTTree->Branch("zhh_phif", &m_zhh_phif, "zhh_phif/F");
  m_pTTree->Branch("zhh_phih", &m_zhh_phih, "zhh_phih/F");
  
  m_pTTree->Branch("zhh_costheta" , &m_zhh_costheta , "zhh_costheta/F");
  m_pTTree->Branch("zhh_costhetaf", &m_zhh_costhetaf, "zhh_costhetaf/F");
  m_pTTree->Branch("zhh_costhetah", &m_zhh_costhetah, "zhh_costhetah/F");

  // 2.b ZHH input
  if (m_saveInputKinematics) {
    m_pTTree->Branch("zhh_l1_e" , &m_zhh_l1_E , "zhh_l1_e/F");
    m_pTTree->Branch("zhh_l1_px", &m_zhh_l1_px, "zhh_l1_px/F");
    m_pTTree->Branch("zhh_l1_py", &m_zhh_l1_py, "zhh_l1_py/F");
    m_pTTree->Branch("zhh_l1_pz", &m_zhh_l1_pz, "zhh_l1_pz/F");

    m_pTTree->Branch("zhh_l2_e" , &m_zhh_l2_E , "zhh_l2_e/F");
    m_pTTree->Branch("zhh_l2_px", &m_zhh_l2_px, "zhh_l2_px/F");
    m_pTTree->Branch("zhh_l2_py", &m_zhh_l2_py, "zhh_l2_py/F");
    m_pTTree->Branch("zhh_l2_pz", &m_zhh_l2_pz, "zhh_l2_pz/F");

    m_pTTree->Branch("zhh_h1_e" , &m_zhh_h1_E , "zhh_h1_e/F");
    m_pTTree->Branch("zhh_h1_px", &m_zhh_h1_px, "zhh_h1_px/F");
    m_pTTree->Branch("zhh_h1_py", &m_zhh_h1_py, "zhh_h1_py/F");
    m_pTTree->Branch("zhh_h1_pz", &m_zhh_h1_pz, "zhh_h1_pz/F");

    m_pTTree->Branch("zhh_h2_e" , &m_zhh_h2_E , "zhh_h2_e/F");
    m_pTTree->Branch("zhh_h2_px", &m_zhh_h2_px, "zhh_h2_px/F");
    m_pTTree->Branch("zhh_h2_py", &m_zhh_h2_py, "zhh_h2_py/F");
    m_pTTree->Branch("zhh_h2_pz", &m_zhh_h2_pz, "zhh_h2_pz/F");
  }

  // 3.a ZZH output
  m_pTTree->Branch("zzh_sigma"  , &m_zzh_sigma  , "zzh_sigma/F");

  m_pTTree->Branch("zzh_sigmalll", &m_zzh_sigmalll, "zzh_sigmalll/F");
  m_pTTree->Branch("zzh_sigmallr", &m_zzh_sigmallr, "zzh_sigmallr/F");

  m_pTTree->Branch("zzh_sigmalrl", &m_zzh_sigmalrl, "zzh_sigmalrl/F");
  m_pTTree->Branch("zzh_sigmalrr", &m_zzh_sigmalrr, "zzh_sigmalrr/F");

  m_pTTree->Branch("zzh_sigmarrr", &m_zzh_sigmarrr, "zzh_sigmarrr/F");
  m_pTTree->Branch("zzh_sigmarrl", &m_zzh_sigmarrl, "zzh_sigmarrl/F");

  m_pTTree->Branch("zzh_sigmarlr", &m_zzh_sigmarlr, "zzh_sigmarlr/F");
  m_pTTree->Branch("zzh_sigmarll", &m_zzh_sigmarll, "zzh_sigmarll/F");

  if (m_zzh_no_z_decay == 1) {
    m_pTTree->Branch("zzh_sigmallz", &m_zzh_sigmallz, "zzh_sigmallz/F");
    m_pTTree->Branch("zzh_sigmalrz", &m_zzh_sigmalrz, "zzh_sigmalrz/F");
    m_pTTree->Branch("zzh_sigmarrz", &m_zzh_sigmarrz, "zzh_sigmarrz/F");
    m_pTTree->Branch("zzh_sigmarlz", &m_zzh_sigmarlz, "zzh_sigmarlz/F");
  }

  m_pTTree->Branch("zzh_mz1"  , &m_zzh_mz1  , "zzh_mz1/F");
  m_pTTree->Branch("zzh_mz2"  , &m_zzh_mz2  , "zzh_mz2/F");
  m_pTTree->Branch("zzh_mzz" , &m_zzh_mzz , "zzh_mzz/F");
  m_pTTree->Branch("zzh_mzzh", &m_zzh_mzzh, "zzh_mzzh/F");
  m_pTTree->Branch("zzh_mh", &m_zzh_mh, "zzh_mh/F");

  m_pTTree->Branch("zzh_phi"   , &m_zzh_phi   , "zzh_phi/F");
  m_pTTree->Branch("zzh_phiz"  , &m_zzh_phiz  , "zzh_phiz/F");
  m_pTTree->Branch("zzh_phiz1f", &m_zzh_phiz1f, "zzh_phiz1f/F");
  m_pTTree->Branch("zzh_phiz2f", &m_zzh_phiz2f, "zzh_phiz2f/F");
  
  m_pTTree->Branch("zzh_costheta"   , &m_zzh_costheta   , "zzh_costheta/F");
  m_pTTree->Branch("zzh_costhetaz"  , &m_zzh_costhetaz  , "zzh_costhetaz/F");
  m_pTTree->Branch("zzh_costhetaz1f", &m_zzh_costhetaz1f, "zzh_costhetaz1f/F");
  m_pTTree->Branch("zzh_costhetaz2f", &m_zzh_costhetaz2f, "zzh_costhetaz2f/F");

  // 3.b ZZH input
  if (m_saveInputKinematics) {
    m_pTTree->Branch("zzh_l1_e" , &m_zzh_l1_E , "zzh_l1_e/F");
    m_pTTree->Branch("zzh_l1_px", &m_zzh_l1_px, "zzh_l1_px/F");
    m_pTTree->Branch("zzh_l1_py", &m_zzh_l1_py, "zzh_l1_py/F");
    m_pTTree->Branch("zzh_l1_pz", &m_zzh_l1_pz, "zzh_l1_pz/F");

    m_pTTree->Branch("zzh_l2_e" , &m_zzh_l2_E , "zzh_l2_e/F");
    m_pTTree->Branch("zzh_l2_px", &m_zzh_l2_px, "zzh_l2_px/F");
    m_pTTree->Branch("zzh_l2_py", &m_zzh_l2_py, "zzh_l2_py/F");
    m_pTTree->Branch("zzh_l2_pz", &m_zzh_l2_pz, "zzh_l2_pz/F");

    m_pTTree->Branch("zzh_z2f1_e" , &m_zzh_z2f1_E , "zzh_z2f1_e/F");
    m_pTTree->Branch("zzh_z2f1_px", &m_zzh_z2f1_px, "zzh_z2f1_px/F");
    m_pTTree->Branch("zzh_z2f1_py", &m_zzh_z2f1_py, "zzh_z2f1_py/F");
    m_pTTree->Branch("zzh_z2f1_pz", &m_zzh_z2f1_pz, "zzh_z2f1_pz/F");

    m_pTTree->Branch("zzh_z2f2_e" , &m_zzh_z2f2_E , "zzh_z2f2_e/F");
    m_pTTree->Branch("zzh_z2f2_px", &m_zzh_z2f2_px, "zzh_z2f2_px/F");
    m_pTTree->Branch("zzh_z2f2_py", &m_zzh_z2f2_py, "zzh_z2f2_py/F");
    m_pTTree->Branch("zzh_z2f2_pz", &m_zzh_z2f2_pz, "zzh_z2f2_pz/F");

    m_pTTree->Branch("zzh_h_e" , &m_zzh_h_E , "zzh_h_e/F");
    m_pTTree->Branch("zzh_h_px", &m_zzh_h_px, "zzh_h_px/F");
    m_pTTree->Branch("zzh_h_py", &m_zzh_h_py, "zzh_h_py/F");
    m_pTTree->Branch("zzh_h_pz", &m_zzh_h_pz, "zzh_h_pz/F");
  }

  // Core configuration
  m_z1_decay_mode = getZDecayModeFromPDG(m_z1_decay_pdg);

  cerr << endl;
  cerr << "ZHH MEM processor initializing with:\n";
  cerr << " H_mass: "<< m_Hmass << "\n";
  cerr << " Z1DecayMode: "<< m_z1_decay_mode << "\n";
  cerr << " Pol_e: "<< pol_e << "\n";
  cerr << " Pol_p: "<< pol_p << "\n";
  cerr << " m_outputTree: " << m_outputTree << "\n";

  _zhh = new LCMEZHH("LCMEZHH", "ZHH", m_Hmass, pol_e, pol_p);
  _zhh->SetZDecayMode(m_z1_decay_mode); // 5 (internal mapping) -> (13) PDG, muon
  _zhh->SetNoHDecay(2);
  _zhh->SetPropagator(1);

  _zzh = new LCMEZZH("LCMEZZH", "ZZH", m_Hmass, pol_e, pol_p);
  _zzh->SetNoZDecay(m_zzh_no_z_decay); // default for ZHH (i.e. no Z decay)
  _zzh->SetNoHDecay(1);
  _zzh->SetPropagator(1);

  // Get matrix element, not diff cross section 
  if (m_mode_me == 1) {
    _zhh->SetMEType(2);
    _zzh->SetMEType(2);
  }
  // zzh setZDecayMode adjusted every run

  streamlog_out(DEBUG) << "   init finished  " << std::endl;
}

void CompareMEProcessor::Clear() 
{
  streamlog_out(DEBUG) << "   Clear called  " << std::endl;

  // Flags to decide whether to save entries for ZHH or ZZH (otherwise they'll all be empty)
  // Important for processes that are not supported by the ME processors (most important, Z->2 boson decay)
  m_zhh_is_set = 0;
  m_zzh_is_set = 0;
  m_error_code = 0;

  // 1. True
  m_true_h1_decay_pdg   = 0;
  m_true_h2_decay_pdg   = 0;
  m_true_z2_decay_pdg   = 0;
  m_passed_preselection = 0;
  m_is_zhh              = 0;
  m_is_zzh              = 0;

  // 2. Event data
  m_h1_decay_pdg  = 0;
  m_h2_decay_pdg  = 0;
  m_z2_decay_pdg  = 0;
  m_z2_decay_mode = 0;

  // 2.a ZHH output
  m_zhh_sigma     = 0.;
  m_zhh_sigmall   = 0.;
  m_zhh_sigmalr   = 0.;
  m_zhh_sigmarl   = 0.;
  m_zhh_sigmarr   = 0.;

  m_zhh_mz        = 0.;
  m_zhh_mhh       = 0.;
  m_zhh_mzhh      = 0.;

  m_zhh_phi       = 0.;
  m_zhh_phif      = 0.;
  m_zhh_phih      = 0.;
  m_zhh_costheta  = 0.;
  m_zhh_costhetaf = 0.;
  m_zhh_costhetah = 0.;

  // 2.b ZHH input
  m_zhh_q2_h1 = 0.;
  m_zhh_q2_h2 = 0.;
  m_zhh_q2_z  = 0.;

  m_zhh_l1_E  = 0.;
  m_zhh_l1_px = 0.;
  m_zhh_l1_py = 0.;
  m_zhh_l1_pz = 0.;

  m_zhh_l2_E  = 0.;
  m_zhh_l2_px = 0.;
  m_zhh_l2_py = 0.;
  m_zhh_l2_pz = 0.;

  m_zhh_h1_E  = 0.;
  m_zhh_h1_px = 0.;
  m_zhh_h1_py = 0.;
  m_zhh_h1_pz = 0.;

  m_zhh_h2_E  = 0.;
  m_zhh_h2_px = 0.;
  m_zhh_h2_py = 0.;
  m_zhh_h2_pz = 0.;


  // 3.a ZZH output
  m_zzh_sigma      = 0.;
  m_zzh_sigmalll   = 0.;
  m_zzh_sigmallz   = 0.;
  m_zzh_sigmallr   = 0.;

  m_zzh_sigmalrl   = 0.;
  m_zzh_sigmalrz   = 0.;
  m_zzh_sigmalrr   = 0.;

  m_zzh_sigmarrr   = 0.;
  m_zzh_sigmarrz   = 0.;
  m_zzh_sigmarrl   = 0.;

  m_zzh_sigmarlr   = 0.;
  m_zzh_sigmarlz   = 0.;
  m_zzh_sigmarll   = 0.;

  m_zzh_mz1       = 0.;
  m_zzh_mz2       = 0.;
  m_zzh_mzz       = 0.;
  m_zzh_mzzh      = 0.;
  
  m_zzh_phi         = 0.;
  m_zzh_phiz        = 0.;
  m_zzh_phiz1f      = 0.;
  m_zzh_phiz2f      = 0.;

  m_zzh_costheta    = 0.;
  m_zzh_costhetaz   = 0.;
  m_zzh_costhetaz1f = 0.;
  m_zzh_costhetaz2f = 0.;

  // 3.b ZZH input
  m_zzh_q2_h  = 0.;
  m_zzh_q2_z1 = 0.;
  m_zzh_q2_z2 = 0.;

  m_zzh_l1_E  = 0.;
  m_zzh_l1_px = 0.;
  m_zzh_l1_py = 0.;
  m_zzh_l1_pz = 0.;

  m_zzh_l2_E  = 0.;
  m_zzh_l2_px = 0.;
  m_zzh_l2_py = 0.;
  m_zzh_l2_pz = 0.;

  m_zzh_z2f1_E  = 0.;
  m_zzh_z2f1_px = 0.;
  m_zzh_z2f1_py = 0.;
  m_zzh_z2f1_pz = 0.;

  m_zzh_z2f2_E  = 0.;
  m_zzh_z2f2_px = 0.;
  m_zzh_z2f2_py = 0.;
  m_zzh_z2f2_pz = 0.;

  m_zzh_h_E  = 0.;
  m_zzh_h_px = 0.;
  m_zzh_h_py = 0.;
  m_zzh_h_pz = 0.;
}

void CompareMEProcessor::processRunHeader( LCRunHeader*  /*run*/) { 
  m_nRun++ ;
}

void CompareMEProcessor::processEvent( EVENT::LCEvent *pLCEvent )
{
  this->Clear();
  
  m_nRun = pLCEvent->getRunNumber();
  m_nEvt = pLCEvent->getEventNumber();
  streamlog_out(DEBUG) << "processing event: " << pLCEvent->getEventNumber() << "  in run: " << pLCEvent->getRunNumber() << endl;

  try {
    // Fetching collections
    LCCollection *inputMCTrueCollection{};
    LCCollection *preselectioncol{};

    streamlog_out(DEBUG) << "        getting true MC collection: " << m_inputMCTrueCollection << std::endl ;
    inputMCTrueCollection = pLCEvent->getCollection( m_inputMCTrueCollection );

    streamlog_out(DEBUG) << "        getting preselection_passed collection: " << m_inputPreSelectionCollection << std::endl ;
    preselectioncol = pLCEvent->getCollection( m_inputPreSelectionCollection );

    // Check for preselection
    m_passed_preselection = preselectioncol->parameters().getIntVal("isPassed");
    if (m_require_presel_pass && !m_passed_preselection)
      return save_evt_with_error_code(ERRORS::PRESELECTION_FAILED_BUT_REQUIRED);

    // Check whether true ZHH or ZZH event
    mcPart_H_if_zhh = (MCParticle*) inputMCTrueCollection->getElementAt(10);
    mcPart_H_if_zzh = (MCParticle*) inputMCTrueCollection->getElementAt(12);

    m_is_zhh = (mcPart_H_if_zhh->getPDG() == 25) && (mcPart_H_if_zzh->getPDG() != 25);
    m_is_zzh = (mcPart_H_if_zzh->getPDG() == 25) && (mcPart_H_if_zhh->getPDG() != 25);

    // Lorentz vectors of final states
    TLorentzVector l1_lortz;
    TLorentzVector l2_lortz;

    TLorentzVector zhh_h1_lortz;
    TLorentzVector zhh_h2_lortz;
    
    TLorentzVector zzh_z2f1_lortz; // only used if m_zzh_no_z_decay == 0
    TLorentzVector zzh_z2f2_lortz; // only used if m_zzh_no_z_decay == 0
    TLorentzVector zzh_z2_lortz;   // only used if m_zzh_no_z_decay == 1
    TLorentzVector zzh_h_lortz;

    // Save truth info about H/Z-decays
    if (m_is_zhh) {
      zhh_h1_decay1 = (MCParticle*) inputMCTrueCollection->getElementAt(12);
      zhh_h2_decay1 = (MCParticle*) inputMCTrueCollection->getElementAt(14);

      m_true_h1_decay_pdg = abs(zhh_h1_decay1->getPDG());
      m_true_h2_decay_pdg = abs(zhh_h2_decay1->getPDG());
    } else {//if (m_is_zzh) {
      zzh_z2f1      = (MCParticle*) inputMCTrueCollection->getElementAt(10);
      zzh_h1_decay1 = (MCParticle*) inputMCTrueCollection->getElementAt(13);

      m_true_z2_decay_pdg = abs(zzh_z2f1->getPDG());
      m_true_h1_decay_pdg = abs(zzh_h1_decay1->getPDG());
    }

    if (m_mode == 0) {
      // MCTruth mode
      m_h1_decay_pdg = m_true_h1_decay_pdg;
      m_h2_decay_pdg = m_true_h2_decay_pdg;
      m_z2_decay_pdg = m_true_z2_decay_pdg;

      // Fetch data from collection holding MCParticle 
      // Get particles of final state
      // Same IDs for final Z1 leptons in both true ZZH and ZHH processes
      l1 = (MCParticle*) inputMCTrueCollection->getElementAt(8);
      l2 = (MCParticle*) inputMCTrueCollection->getElementAt(9);

      l1_lortz = v4(l1);
      l2_lortz = v4(l2);

      if (m_is_zhh) {
        // Assumming ZHH
        zhh_h1 = (MCParticle*) inputMCTrueCollection->getElementAt(10);
        zhh_h2 = (MCParticle*) inputMCTrueCollection->getElementAt(11);

        zhh_h1_lortz = v4(zhh_h1);
        zhh_h2_lortz = v4(zhh_h2);

        m_zhh_is_set = 1;

        // Assuming ZZH
        if (m_zzh_no_z_decay == 0) {
          // Pretend that decay products of H1 are decay products of Z2 in ZZH
          zhh_h1_decay2 = (MCParticle*) inputMCTrueCollection->getElementAt(13);
          
          zzh_z2f1_lortz = v4(zhh_h1_decay1);
          zzh_z2f2_lortz = v4(zhh_h1_decay2);

          m_z2_decay_mode = getZDecayModeFromPDG(m_h1_decay_pdg);

          if (m_z2_decay_mode > 0) {
            zzh_h_lortz = zhh_h2_lortz; // identify H2 in ZHH as H in ZZH
            m_zzh_is_set = 1;
          }
        } else if (m_zzh_no_z_decay == 1) {
          zzh_z2_lortz = zhh_h1_lortz;
          zzh_h_lortz  = zhh_h2_lortz;
          m_zzh_is_set = 1;
        }
        
      } else {//if (m_is_zzh) {
        // Assumming ZZH
        //          zzh_z2f1 already assigned above, in truth info block
        zzh_z2f2 = (MCParticle*) inputMCTrueCollection->getElementAt(11);
        zzh_h    = (MCParticle*) inputMCTrueCollection->getElementAt(12);

        zzh_z2f1_lortz = v4(zzh_z2f1->getPDG() > 0 ? zzh_z2f1 : zzh_z2f2);
        zzh_z2f2_lortz = v4(zzh_z2f1->getPDG() > 0 ? zzh_z2f2 : zzh_z2f1);

        zzh_h_lortz    = v4(zzh_h);

        if (m_zzh_no_z_decay == 0) {
          m_z2_decay_mode = getZDecayModeFromPDG(m_z2_decay_pdg);
        } else if (m_zzh_no_z_decay == 1) {
          zzh_z2_lortz = zzh_z2f1_lortz + zzh_z2f2_lortz;
        }

        m_zzh_is_set = 1;

        // Assuming ZHH
        zhh_h1_lortz = zzh_z2f1_lortz + zzh_z2f2_lortz;
        zhh_h2_lortz = zzh_h_lortz;
        
        m_zhh_is_set = 1;
      }
    } else if (m_mode == 1 || m_mode == 2) {
      // Jet reconstruction mode
      // 1 for RefinedJets (or any m_inputJetCollection) with inputHiggsPair; or 2 for TrueJet (requires output as from Misclustering) in https://github.com/nVentis/ZHH
      // As of now, assumes inputJetCol, inputLepPair, inputHiggsPair and inputHdecayMode

      // To be filled in both modes
      vector<ReconstructedParticle*> reco_jets; // not using "jets" because its reserved in TrueJet_Parser
      vector<int> perm; 

      // TrueJet helper
      this->getall( pLCEvent );
      DelMe delme(std::bind(&CompareMEProcessor::delall, this));

      // Fetch hadronic jets
      if (m_mode == 2) {
        // Check if JetMatching collection exists
        const vector<string> *coll_names = pLCEvent->getCollectionNames();
        if (std::find(coll_names->begin(), coll_names->end(), m_inputJetMatchingCollection) == coll_names->end())
          return save_evt_with_error_code(ERRORS::NO_JET_MATCHING_COLLECTION);

        // Fetching collections
        streamlog_out(DEBUG) << " getting JetMatching collection: " << m_inputJetMatchingCollection << std::endl ;
        const EVENT::LCParameters& jm_params = pLCEvent->getCollection( m_inputJetMatchingCollection )->getParameters();

        vector<int> temp;
        if (jm_params.getNInt(std::string("h1tjet1id")) > 0) {
          // Supported in newer versions
          temp.push_back(jm_params.getIntVal("h1tjet1id"));
          temp.push_back(jm_params.getIntVal("h1tjet2id"));
          temp.push_back(jm_params.getIntVal("h2tjet1id"));
          temp.push_back(jm_params.getIntVal("h2tjet2id"));

          for (int j = 0; j < 4; j++) {
            ReconstructedParticle* true_jet = (ReconstructedParticle*) jet(temp[j]);

            reco_jets.push_back(true_jet);
            perm.push_back(j); // jets are already matched by true_jet_index
          }
        } else {
          temp.push_back(jm_params.getIntVal("h1jet1id"));
          temp.push_back(jm_params.getIntVal("h1jet2id"));
          temp.push_back(jm_params.getIntVal("h2jet1id"));
          temp.push_back(jm_params.getIntVal("h2jet2id"));

          for (int j = 0; j < 4; j++) {
            char key[12];
            sprintf(key, "reco2true%d", temp[j]);

            int true_jet_index = jm_params.getIntVal(key);

            ReconstructedParticle* true_jet = (ReconstructedParticle*) jet(true_jet_index);
            reco_jets.push_back(true_jet);
            perm.push_back(j); // jets are already matched by true_jet_index
          }
        }
      } else if (m_mode == 1) {
        // Fetching collections
        streamlog_out(DEBUG) << " getting jet collection: " << m_inputJetCollection << std::endl;
        LCCollection *inputJetCol = pLCEvent->getCollection( m_inputJetCollection );

        streamlog_out(DEBUG) << " getting higgs pair: " << m_inputHiggsPairCollection << std::endl;
        LCCollection *inputHiggsPair = pLCEvent->getCollection( m_inputHiggsPairCollection );

        // Get Higgs decay jets
        const EVENT::LCParameters& higgsParams = inputHiggsPair->getParameters();

        // Fetch jets and isolated leptons
        for (int i=0; i < 4; ++i) {
          ReconstructedParticle* reco_jet = (ReconstructedParticle*) inputJetCol->getElementAt(i);
          reco_jets.push_back(reco_jet);
        }

        // Check if jet pairing parameters exist in higgs pair; otherwise try permutations to check which pairing was used  
        if (higgsParams.getNInt(std::string("h1jet1id")) == 1) {
          // Jet pairing saved in collection (newer version)
          perm.push_back(higgsParams.getIntVal("h1jet1id"));
          perm.push_back(higgsParams.getIntVal("h1jet2id"));
          perm.push_back(higgsParams.getIntVal("h2jet1id"));
          perm.push_back(higgsParams.getIntVal("h2jet2id"));
        } else {
          vector<vector<int>> perms {
            {0,1,2,3},
            {0,2,1,3},
            {0,3,1,2}
          };

          float chi2min = 99999. ;

          // Pair first jet with three others and check where match is best
          int best_idx = 0;
          for (unsigned int i = 0; i < perms.size(); i++) {
            float m1 = inv_mass(reco_jets[perms[i][0]], reco_jets[perms[i][1]]);
            float chi2 = (m1-125)*(m1-125);

            if (chi2 < chi2min) {
              chi2min = chi2;
              best_idx = i;
            }
          }

          perm = perms[best_idx];
        }
      }

      // Fetch leptons (errors 101 and 201)
      if (m_lepton_mode == 2 || (m_lepton_mode == 0 && m_mode == 2)) {
        unsigned int n_matching_l_jets = 0;
        unsigned int n_matching_lb_jets = 0;
        for (int i_jet = 0 ; i_jet < this->njets() ; i_jet++ ) {
          if ( type_jet( i_jet ) == 2 ) { // 2 = leptonic jets
            int jet_pdg = jet(i_jet)->getParticleIDs()[0]->getPDG();
            if ( jet_pdg == m_z1_decay_pdg ) {
              l1 = (ReconstructedParticle*) jet(i_jet);
              n_matching_l_jets++;
            } else if ( jet_pdg == -m_z1_decay_pdg) {
              l2 = (ReconstructedParticle*) jet(i_jet);
              n_matching_lb_jets++;
            }
          }
        }

        if (n_matching_l_jets != n_matching_lb_jets || n_matching_lb_jets != 1)
          return save_evt_with_error_code(ERRORS::INCOMPLETE_TRUEJET_LEPTON_PAIR);

      } else if (m_lepton_mode == 1 || (m_lepton_mode == 0 && m_mode == 1)) {
        streamlog_out(DEBUG) << " getting lepton pair: " << m_inputLepPairCollection << std::endl;
        LCCollection *inputLepPair = pLCEvent->getCollection( m_inputLepPairCollection );

        if (inputLepPair->getNumberOfElements() != 2)
          return save_evt_with_error_code(ERRORS::INCOMPLETE_RECO_LEPTON_PAIR);

        l1 = (ReconstructedParticle*) inputLepPair->getElementAt(0);
        l2 = (ReconstructedParticle*) inputLepPair->getElementAt(1)

        l1_lortz = v4(l1);
        l2_lortz = v4(l2);

        free(l1);
        free(l2);
      }

      // CHEATED: Extract information about decays
      int both_to_b = 0;
      int both_to_c = 0;
      if (m_is_zhh) {
        both_to_b = (m_true_h1_decay_pdg == m_true_h2_decay_pdg) && (m_true_h1_decay_pdg == 5);
        both_to_c = (m_true_h1_decay_pdg == m_true_h2_decay_pdg) && (m_true_h1_decay_pdg == 4);
      } else if (m_is_zzh) {
        both_to_b = (m_true_h1_decay_pdg == m_true_z2_decay_pdg) && (m_true_h1_decay_pdg == 5);
        both_to_c = (m_true_h1_decay_pdg == m_true_z2_decay_pdg) && (m_true_h1_decay_pdg == 4);
      }
      
      // Necessary for ME calculation: b/c jets (only, because we cannot tell apart which [reco] jet might be which?) AND 2 isolated leptons (i.e. one LeptonPair)
      // cerr << " both_to_b " << both_to_b << " both_to_c " << both_to_c << " lep_pair " << inputLepPair->getNumberOfElements() << endl;
      if (!(both_to_b || both_to_c))
        return save_evt_with_error_code(ERRORS::NEITHER_BBBB_NOR_CCCC);
    
      streamlog_out(MESSAGE) << "processEvent : using Higgs jet pairing " << perm[0] << perm[1] << perm[2] << perm[3] << std::endl;

      // Assuming ZHH
      zzh_z2f1_lortz = v4(reco_jets[perm[2]]);
      zzh_z2f2_lortz = v4(reco_jets[perm[3]]);
      zzh_h_lortz    = v4(reco_jets[perm[0]]) + v4(reco_jets[perm[1]]);

      m_z2_decay_pdg  = both_to_b ? 5 : 4; // PDGs: bottom->5, charm->4
      m_h1_decay_pdg  =  both_to_b ? 5 : 4;
      m_z2_decay_mode = getZDecayModeFromPDG(m_z2_decay_pdg);

      m_zhh_is_set = 1;

      // Assuming ZHH
      zhh_h1_lortz = zzh_z2f1_lortz + zzh_z2f2_lortz;
      zhh_h2_lortz = zzh_h_lortz;

      m_h2_decay_pdg = both_to_b ? 5 : 4;
      
      m_zzh_is_set = 1;
    }

    // ZHH
    if (m_zhh_is_set) {
      // Set final states
      TLorentzVector zhh_lortz[4] = {l1_lortz, l2_lortz, zhh_h1_lortz, zhh_h2_lortz};

      _zhh->SetMomentumFinal(zhh_lortz);
      // Zdecay unchanged

      // Output
      Int_t vHelLL[2] = {-1,-1};
      Int_t vHelLR[2] = {-1, 1};
      Int_t vHelRL[2] = { 1,-1};
      Int_t vHelRR[2] = { 1, 1};
      
      m_zhh_mz   = TMath::Sqrt(_zhh->GetQ2Z());
      m_zhh_mhh  = TMath::Sqrt(_zhh->GetQ2HH());
      m_zhh_mzhh = TMath::Sqrt(_zhh->GetQ2ZHH());

      m_zhh_sigma   = _zhh->GetMatrixElement2();
      m_zhh_sigmall = _zhh->GetMatrixElement2(vHelLL);
      m_zhh_sigmalr = _zhh->GetMatrixElement2(vHelLR);
      m_zhh_sigmarl = _zhh->GetMatrixElement2(vHelRL);
      m_zhh_sigmarr = _zhh->GetMatrixElement2(vHelRR);

      m_zhh_phi       = _zhh->GetPhi();
      m_zhh_phif      = _zhh->GetPhiF();
      m_zhh_phih      = _zhh->GetPhiH();
      m_zhh_costheta  = _zhh->GetCosTheta();
      m_zhh_costhetaf = _zhh->GetCosThetaF();
      m_zhh_costhetah = _zhh->GetCosThetaH();

      // Input
      m_zhh_q2_h1 = _zhh->GetQ2H1();
      m_zhh_q2_h2 = _zhh->GetQ2H2();
      m_zhh_q2_z  = _zhh->GetQ2Z();

      m_zhh_h1_E  = zhh_h1_lortz.E();
      m_zhh_h1_px = zhh_h1_lortz.Px();
      m_zhh_h1_py = zhh_h1_lortz.Py();
      m_zhh_h1_pz = zhh_h1_lortz.Pz();

      m_zhh_h2_E  = zhh_h2_lortz.E();
      m_zhh_h2_px = zhh_h2_lortz.Px();
      m_zhh_h2_py = zhh_h2_lortz.Py();
      m_zhh_h2_pz = zhh_h2_lortz.Pz();
    }

    // ZZH
    if (m_zzh_is_set) {
      // Possible helicity initial and final states
      Int_t vHelLLL[3] = {-1,-1,-1};
      Int_t vHelLLR[3] = {-1,-1, 1};
      Int_t vHelLRL[3] = {-1, 1,-1};
      Int_t vHelLRR[3] = {-1, 1, 1};
      Int_t vHelRLL[3] = { 1,-1,-1};
      Int_t vHelRLR[3] = { 1,-1, 1};
      Int_t vHelRRL[3] = { 1, 1,-1};
      Int_t vHelRRR[3] = { 1, 1, 1};
      
      // Set final states
      if (m_zzh_no_z_decay == 0) {
        TLorentzVector zzh_lortz[5] = {l1_lortz, l2_lortz, zzh_z2f1_lortz, zzh_z2f2_lortz, zzh_h_lortz};
        
        _zzh->SetZDecayMode(m_z1_decay_mode, m_z2_decay_mode);
        _zzh->SetMomentumFinal(zzh_lortz);
      } else if (m_zzh_no_z_decay == 1) {
        TLorentzVector zzh_lortz[4] = {l1_lortz, l2_lortz, zzh_z2_lortz, zzh_h_lortz};

        _zzh->SetZDecayMode(m_z1_decay_mode);
        _zzh->SetMomentumFinalNoZdecay(zzh_lortz);

        Int_t vHelLLZ[3] = {-1,-1, 0};
        Int_t vHelLRZ[3] = {-1, 1, 0};
        Int_t vHelRLZ[3] = { 1,-1, 0};
        Int_t vHelRRZ[3] = { 1, 1, 0};

        m_zzh_sigmallz = _zzh->GetMatrixElement2(vHelLLZ);
        m_zzh_sigmalrz = _zzh->GetMatrixElement2(vHelLRZ);
        m_zzh_sigmarlz = _zzh->GetMatrixElement2(vHelRLZ);
        m_zzh_sigmarrz = _zzh->GetMatrixElement2(vHelRRZ);
      }

      m_zzh_mz1  = TMath::Sqrt(_zzh->GetQ2Z1());
      m_zzh_mz2  = TMath::Sqrt(_zzh->GetQ2Z2());
      m_zzh_mzz  = TMath::Sqrt(_zzh->GetQ2ZZ());
      m_zzh_mzzh = TMath::Sqrt(_zzh->GetQ2ZZH());
      m_zzh_mh   = TMath::Sqrt(_zzh->GetQ2H());

      m_zzh_sigma    = _zzh->GetMatrixElement2();
      m_zzh_sigmalll = _zzh->GetMatrixElement2(vHelLLL);
      m_zzh_sigmallr = _zzh->GetMatrixElement2(vHelLLR);
      m_zzh_sigmalrl = _zzh->GetMatrixElement2(vHelLRL);
      m_zzh_sigmalrr = _zzh->GetMatrixElement2(vHelLRR);
      m_zzh_sigmarll = _zzh->GetMatrixElement2(vHelRLL);
      m_zzh_sigmarlr = _zzh->GetMatrixElement2(vHelRLR);
      m_zzh_sigmarrl = _zzh->GetMatrixElement2(vHelRRL);
      m_zzh_sigmarrr = _zzh->GetMatrixElement2(vHelRRR);

      m_zzh_phi    = _zzh->GetPhi();
      m_zzh_phiz1f = _zzh->GetPhiZ1F();
      m_zzh_phiz2f = _zzh->GetPhiZ2F();
      m_zzh_phiz   = _zzh->GetPhiZ();
      m_zzh_costheta    = _zzh->GetCosTheta();
      m_zzh_costhetaz1f = _zzh->GetCosThetaZ1F();
      m_zzh_costhetaz2f = _zzh->GetCosThetaZ2F();
      m_zzh_costhetaz   = _zzh->GetCosThetaZ();

      // Input
      m_zzh_q2_z1 = _zzh->GetQ2Z1();
      m_zzh_q2_z2 = _zzh->GetQ2Z2();
      m_zzh_q2_h  = _zzh->GetQ2H();

      if (m_zzh_no_z_decay == 0) {
        m_zzh_z2f1_E  = zzh_z2f1_lortz.E();
        m_zzh_z2f1_px = zzh_z2f1_lortz.Px();
        m_zzh_z2f1_py = zzh_z2f1_lortz.Py();
        m_zzh_z2f1_pz = zzh_z2f1_lortz.Pz();
      } else if (m_zzh_no_z_decay == 1) {
        m_zzh_z2f1_E  = zzh_z2_lortz.E();
        m_zzh_z2f1_px = zzh_z2_lortz.Px();
        m_zzh_z2f1_py = zzh_z2_lortz.Py();
        m_zzh_z2f1_pz = zzh_z2_lortz.Pz();
      }

      m_zzh_z2f2_E  = zzh_z2f2_lortz.E();
      m_zzh_z2f2_px = zzh_z2f2_lortz.Px();
      m_zzh_z2f2_py = zzh_z2f2_lortz.Py();
      m_zzh_z2f2_pz = zzh_z2f2_lortz.Pz();

      m_zzh_h_E  = zzh_h_lortz.E();
      m_zzh_h_px = zzh_h_lortz.Px();
      m_zzh_h_py = zzh_h_lortz.Py();
      m_zzh_h_pz = zzh_h_lortz.Pz();
    }

    // ZHH and ZZH input (shared l1 and l2 leptons)
    m_zhh_l1_E  = l1_lortz.E();
    m_zhh_l1_px = l1_lortz.Px();
    m_zhh_l1_py = l1_lortz.Py();
    m_zhh_l1_pz = l1_lortz.Pz();

    m_zhh_l2_E  = l2_lortz.E();
    m_zhh_l2_px = l2_lortz.Px();
    m_zhh_l2_py = l2_lortz.Py();
    m_zhh_l2_pz = l2_lortz.Pz();

    m_zzh_l1_E  = l1_lortz.E();
    m_zzh_l1_px = l1_lortz.Px();
    m_zzh_l1_py = l1_lortz.Py();
    m_zzh_l1_pz = l1_lortz.Pz();

    m_zzh_l2_E  = l2_lortz.E();
    m_zzh_l2_px = l2_lortz.Px();
    m_zzh_l2_py = l2_lortz.Py();
    m_zzh_l2_pz = l2_lortz.Pz();

    m_pTTree->Fill();
    
  } catch(DataNotAvailableException &e) {
    streamlog_out(MESSAGE) << "processEvent : Input collections not found in event " << m_nEvt << std::endl;
  }

}

void CompareMEProcessor::check()
{
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}

void CompareMEProcessor::save_evt_with_error_code(int error_code)
{
  m_error_code = error_code;
  m_pTTree->Fill();
}


void CompareMEProcessor::end()
{
  m_pTFile->cd();
  m_pTTree->Write();
  m_pTFile->Close();
  delete m_pTFile;
}

int CompareMEProcessor::getZDecayModeFromPDG(int pdg)
{
  Int_t pdgs[11] = {12, 14, 16, 11, 13, 15, 2, 4, 1, 3 ,  5};
  Int_t zdms[11] = {1 , 2,  3,  4,  5,  6,  7, 8, 9, 10, 11};

  for (unsigned int i = 0; i < 11; i++) {
    if (pdgs[i] == pdg) {
      return zdms[i];
    }
  }

  return 0;
}