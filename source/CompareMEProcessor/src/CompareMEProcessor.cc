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
#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>
#include <IMPL/ReconstructedParticleImpl.h>
#include <UTIL/PIDHandler.h>

using namespace lcio ;
using namespace marlin ;
using namespace lcme;

template<class T>
TLorentzVector v4(T* p){
  return TLorentzVector( p->getMomentum(), p->getEnergy() );
}

/**
 * Wrapper for TrueJet to TLorentzVector
 * For TrueJet p4, p4[0] is energy, unlike for ROOT
*/
TLorentzVector v4(const Double_t * x0){
  return TLorentzVector( x0[1], x0[2], x0[3], x0[0] );
}

template<class T>
void set_v4(TLorentzVector &t, T* p){
  t.SetPx(p->getMomentum()[0]);
  t.SetPy(p->getMomentum()[1]);
  t.SetPz(p->getMomentum()[2]);
  t.SetE(p->getEnergy());
}

template<class T>
double inv_mass(T* p1, T* p2){
  double e = p1->getEnergy()+p2->getEnergy() ;
  double px = p1->getMomentum()[0]+p2->getMomentum()[0];
  double py = p1->getMomentum()[1]+p2->getMomentum()[1];
  double pz = p1->getMomentum()[2]+p2->getMomentum()[2];
  return( sqrt( e*e - px*px - py*py - pz*pz  ) );
}

CompareMEProcessor aCompareMEProcessor ;

CompareMEProcessor::CompareMEProcessor() :

  Processor("CompareMEProcessor"),
  m_nRun(0),
  m_nEvt(0),
  m_mode_me(1), // 0-> dsigma, 1-> ME^2

  parton1(0., 0., 0., 0.),
  parton2(0., 0., 0., 0.),
  parton3(0., 0., 0., 0.),
  parton4(0., 0., 0., 0.),
  true_lep1(0., 0., 0., 0.),
  true_lep2(0., 0., 0., 0.),

  jet1(0., 0., 0., 0.),
  jet2(0., 0., 0., 0.),
  jet3(0., 0., 0., 0.),
  jet4(0., 0., 0., 0.),
  lep1(0., 0., 0., 0.),
  lep2(0., 0., 0., 0.)
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
				 "JetMatchingSigCollection",
				 "Name of collection c parameters about jet matching for the signal hypothesis",
				 m_inputJetMatchingCollection,
				 std::string("JetMatchingSigDurham")
				 );

  registerProcessorParameter("RequirePreselectionPass",
        "whether (1) or not to skip (1) events that did not pass preselection",
        m_require_presel_pass,
        int(0)
        );

  registerProcessorParameter("RequireNToXPass",
        "how many ICNs are required to decay to ",
        m_require_ntox_pass,
        int(0)
        );

  registerProcessorParameter("RequireNToX",
        "which PDG ICNs are required to decay to; defaults to 5/bottom",
        m_require_ntox,
        int(5)
        );
  /*
  registerProcessorParameter("RequireMCEnergyConservation",
        "whether (1) or not to skip (1) events that did not pass preselection",
        m_require_mc_energycnsv,
        int(0)
        );*/

	registerProcessorParameter("Z1DecayPDG",
        "MEM processor mode of decay of Z1 (=Z for ZHH, Z1 for ZZH)",
        m_z1_decay_pdg,
        int(13)
        );

  registerProcessorParameter("Mode",
        "calculate MEs based on 0: True/MCParticleSkimmed; (following both using JetClusteringSig matching) 1: RefinedJets; 2: TrueJet (with matching of Misclustering processor in ZHH github); 3: TrueJet (with matching according to TrueJet)",
        m_mode,
        int(0)
        );

  registerProcessorParameter("LeptonMode",
        "which data source to assume for the Z->mu+mu- decay; -1: same as Mode, 0:MCTruth, 1: LeptonPair, 2: TrueJet (deprecated)",
        m_lepton_mode,
        int(-1)
        );

  registerProcessorParameter("TrueJetMode",
        "0 for true, 1 for seen",
        m_truejet_mode,
        int(0)
        );

  registerProcessorParameter("SaveTransferEnergies",
        "0: dont save parton, jet and lepton energies; 1: do it",
        m_saveTransferEnergies,
        int(1)
        );

  registerProcessorParameter("SaveTransferKinematics",
        "0: dont save parton, jets and lepton momenta; 1: do it",
        m_saveTransferKinematics,
        int(1)
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

  registerInputCollection( LCIO::LCRELATION,
                            "RecoMCTruthLink",
                            "Name of the RecoMCTruthLink input collection"  ,
                            _recoMCTruthLink,
                            string("RecoMCTruthLink")
                            );
}

void CompareMEProcessor::init()
{
  streamlog_out(DEBUG) << "   init called  " << std::endl;
  this->Clear();

  // Configuration
  Double_t pol_e = 1.; // this->parameters()->getFloatVal("beamPol1");
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

  // For Reco and TrueJet mode, save the region and region_icns parameters of the Misclustering processor
  if (m_mode > 0) {
    m_pTTree->Branch("misclustering_region", &m_misclustering_region, "misclustering_region/I");
    m_pTTree->Branch("misclustering_region_icns", &m_misclustering_region_icns, "misclustering_region_icns/I");

    m_pTTree->Branch("efrac1_reco", &m_efrac1_reco, "efrac1_icn_reco/F");
    m_pTTree->Branch("efrac2_reco", &m_efrac2_reco, "efrac1_icn_reco/F");
    m_pTTree->Branch("efrac1_true", &m_efrac1_true, "efrac1_icn_true/F");
    m_pTTree->Branch("efrac2_true", &m_efrac2_true, "efrac1_icn_true/F");

    m_pTTree->Branch("efrac1_icn_reco", &m_efrac1_icn_reco, "efrac1_icn_reco/F");
    m_pTTree->Branch("efrac2_icn_reco", &m_efrac2_icn_reco, "efrac1_icn_reco/F");
    m_pTTree->Branch("efrac1_icn_true", &m_efrac1_icn_true, "efrac1_icn_true/F");
    m_pTTree->Branch("efrac2_icn_true", &m_efrac2_icn_true, "efrac1_icn_true/F");
  }

  m_pTTree->Branch("parton1_pdg", &m_parton1_pdg);
  m_pTTree->Branch("parton2_pdg", &m_parton2_pdg);
  m_pTTree->Branch("parton3_pdg", &m_parton3_pdg);
  m_pTTree->Branch("parton4_pdg", &m_parton4_pdg);

  if (m_saveTransferEnergies > 0) {
    m_pTTree->Branch("parton1_e", &m_parton1_e);
    m_pTTree->Branch("parton2_e", &m_parton2_e);
    m_pTTree->Branch("parton3_e", &m_parton3_e);
    m_pTTree->Branch("parton4_e", &m_parton4_e);

    m_pTTree->Branch("true_lep1_e", &m_true_lep1_e);
    m_pTTree->Branch("true_lep2_e", &m_true_lep2_e);

    m_pTTree->Branch("jet1_e", &m_jet1_e);
    m_pTTree->Branch("jet2_e", &m_jet2_e);
    m_pTTree->Branch("jet3_e", &m_jet3_e);
    m_pTTree->Branch("jet4_e", &m_jet4_e);

    m_pTTree->Branch("lep1_e", &m_lep1_e);
    m_pTTree->Branch("lep2_e", &m_lep2_e);
  }

  if (m_saveTransferKinematics > 0) {
    m_pTTree->Branch("parton1_px", &m_parton1_px);
    m_pTTree->Branch("parton1_py", &m_parton1_py);
    m_pTTree->Branch("parton1_pz", &m_parton1_pz);

    m_pTTree->Branch("parton2_px", &m_parton2_px);
    m_pTTree->Branch("parton2_py", &m_parton2_py);
    m_pTTree->Branch("parton2_pz", &m_parton2_pz);

    m_pTTree->Branch("parton3_px", &m_parton3_px);
    m_pTTree->Branch("parton3_py", &m_parton3_py);
    m_pTTree->Branch("parton3_pz", &m_parton3_pz);

    m_pTTree->Branch("parton4_px", &m_parton4_px);
    m_pTTree->Branch("parton4_py", &m_parton4_py);
    m_pTTree->Branch("parton4_pz", &m_parton4_pz);

    m_pTTree->Branch("true_lep1_px", &m_true_lep1_px);
    m_pTTree->Branch("true_lep1_py", &m_true_lep1_py);
    m_pTTree->Branch("true_lep1_pz", &m_true_lep1_pz);

    m_pTTree->Branch("true_lep2_px", &m_true_lep2_px);
    m_pTTree->Branch("true_lep2_py", &m_true_lep2_py);
    m_pTTree->Branch("true_lep2_pz", &m_true_lep2_pz);

    m_pTTree->Branch("jet1_px", &m_jet1_px);
    m_pTTree->Branch("jet1_py", &m_jet1_py);
    m_pTTree->Branch("jet1_pz", &m_jet1_pz);

    m_pTTree->Branch("jet2_px", &m_jet2_px);
    m_pTTree->Branch("jet2_py", &m_jet2_py);
    m_pTTree->Branch("jet2_pz", &m_jet2_pz);

    m_pTTree->Branch("jet3_px", &m_jet3_px);
    m_pTTree->Branch("jet3_py", &m_jet3_py);
    m_pTTree->Branch("jet3_pz", &m_jet3_pz);

    m_pTTree->Branch("jet4_px", &m_jet4_px);
    m_pTTree->Branch("jet4_py", &m_jet4_py);
    m_pTTree->Branch("jet4_pz", &m_jet4_pz);

    m_pTTree->Branch("lep1_px", &m_lep1_px);
    m_pTTree->Branch("lep1_py", &m_lep1_py);
    m_pTTree->Branch("lep1_pz", &m_lep1_pz);

    m_pTTree->Branch("lep2_px", &m_lep2_px);
    m_pTTree->Branch("lep2_py", &m_lep2_py);
    m_pTTree->Branch("lep2_pz", &m_lep2_pz);
  }

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

  // Core configuration
  m_z1_decay_mode = getZDecayModeFromPDG(m_z1_decay_pdg);

  std::cerr << std::endl;
  std::cerr << "ZHH MEM processor initializing with:\n";
  std::cerr << " H_mass: "<< m_Hmass << "\n";
  std::cerr << " Z1DecayMode: "<< m_z1_decay_mode << "\n";
  std::cerr << " Pol_e: "<< pol_e << "\n";
  std::cerr << " Pol_p: "<< pol_p << "\n";
  std::cerr << " m_outputTree: " << m_outputTree << "\n";

  _zhh = new LCMEZHH("LCMEZHH", "ZHH", m_Hmass, pol_e, pol_p);
  _zhh->SetZDecayMode(m_z1_decay_mode); // 5 (internal mapping) -> (13) PDG, muon
  _zhh->SetNoHDecay(2);
  _zhh->SetPropagator(1);

  _zzh = new LCMEZZH("LCMEZZH", "ZZH", m_Hmass, pol_e, pol_p);
  _zzh->SetNoZDecay(0);
  _zzh->SetNoHDecay(1);
  _zzh->SetPropagator(1);

  // Get matrix element, not diff cross section 
  if (m_mode_me == 1) {
    _zhh->SetMEType(2);
    _zzh->SetMEType(2);
  }
  // zzh setZDecayMode adjusted every run

  streamlog_out(DEBUG) << "init finished  " << std::endl;
}

void CompareMEProcessor::Clear() 
{
  streamlog_out(DEBUG) << "clear called  " << std::endl;

  // Only delete LCRelationNavigator onstances in TrueJet_Parser if these objects have been instantiated before with "new" 
  if (m_error_code != ERRORS::PRESELECTION_FAILED_BUT_REQUIRED &&
      m_error_code != ERRORS::NTOX_FAILED_BUT_REQUIRED &&
      m_error_code != ERRORS::TRUEJET_NO_JETS &&
      m_error_code != ERRORS::TRUEJET_NULL)
    this->delall();

  // Flags to decide whether to save entries for ZHH or ZZH (otherwise they'll all be empty)
  // Important for processes that are not supported by the ME processors (most important, Z->2 boson decay)
  m_zhh_is_set = 0;
  m_zzh_is_set = 0;
  m_error_code = 0;

  // Energy+impulse transfer
  parton1.SetPxPyPzE(0, 0, 0, 0);
  parton2.SetPxPyPzE(0, 0, 0, 0);
  parton3.SetPxPyPzE(0, 0, 0, 0);
  parton4.SetPxPyPzE(0, 0, 0, 0);

  true_lep1.SetPxPyPzE(0, 0, 0, 0);
  true_lep2.SetPxPyPzE(0, 0, 0, 0);

  jet1.SetPxPyPzE(0, 0, 0, 0);
  jet2.SetPxPyPzE(0, 0, 0, 0);
  jet3.SetPxPyPzE(0, 0, 0, 0);
  jet4.SetPxPyPzE(0, 0, 0, 0);

  lep1.SetPxPyPzE(0, 0, 0, 0);
  lep2.SetPxPyPzE(0, 0, 0, 0);

  m_parton1_pdg = 0;
  m_parton2_pdg = 0;
  m_parton3_pdg = 0;
  m_parton4_pdg = 0;

  m_parton1_e  = m_parton2_e  = m_parton3_e  = m_parton4_e  = 0.;
  m_parton1_px = m_parton2_px = m_parton3_px = m_parton4_px = 0.;
  m_parton1_py = m_parton2_py = m_parton3_py = m_parton4_py = 0.;
  m_parton1_pz = m_parton2_pz = m_parton3_pz = m_parton4_pz = 0.;

  m_true_lep1_e  = m_true_lep2_e  = 0.;
  m_true_lep1_px = m_true_lep2_px = 0.;
  m_true_lep1_py = m_true_lep2_py = 0.;
  m_true_lep1_pz = m_true_lep2_pz = 0.;

  m_jet1_e  = m_jet2_e  = m_jet3_e  = m_jet4_e  = 0.;
  m_jet1_px = m_jet2_px = m_jet3_px = m_jet4_px = 0.;
  m_jet1_py = m_jet2_py = m_jet3_py = m_jet4_py = 0.;
  m_jet1_pz = m_jet2_pz = m_jet3_pz = m_jet4_pz = 0.;

  m_lep1_e  = m_lep2_e  = 0.;
  m_lep1_px = m_lep2_px = 0.;
  m_lep1_py = m_lep2_py = 0.;
  m_lep1_pz = m_lep2_pz = 0.;

  // 1. True
  m_true_h1_decay_pdg         = 0;
  m_true_h2_decay_pdg         = 0;
  m_true_z2_decay_pdg         = 0;
  m_passed_preselection       = 0;
  m_is_zhh                    = 0;
  m_is_zzh                    = 0;
  m_misclustering_region      = -1;
  m_misclustering_region_icns = -1;

  m_efrac1_reco = 0.;
  m_efrac2_reco = 0.;
  m_efrac1_true = 0.;
  m_efrac2_true = 0.;

  m_efrac1_icn_reco = 0.;
  m_efrac2_icn_reco = 0.;
  m_efrac1_icn_true = 0.;
  m_efrac2_icn_true = 0.;

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

  m_zhh_q2_h1     = 0.;
  m_zhh_q2_h2     = 0.;
  m_zhh_q2_z      = 0.;

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

  m_zzh_q2_h        = 0.;
  m_zzh_q2_z1       = 0.;
  m_zzh_q2_z2       = 0.;
}

void CompareMEProcessor::processRunHeader( LCRunHeader*  /*run*/) { 
  m_nRun++ ;
}

// Writes the data in four-momenta into their corresponding float/double values for writing them to ROOT
void CompareMEProcessor::transfer_four_momenta()
{
  m_parton1_e = parton1.E();
  m_parton2_e = parton2.E();
  m_parton3_e = parton3.E();
  m_parton4_e = parton4.E();

  m_parton1_px = parton1.Px();
  m_parton1_py = parton1.Py();
  m_parton1_pz = parton1.Pz();

  m_parton2_px = parton2.Px();
  m_parton2_py = parton2.Py();
  m_parton2_pz = parton2.Pz();

  m_parton3_px = parton3.Px();
  m_parton3_py = parton3.Py();
  m_parton3_pz = parton3.Pz();

  m_parton4_px = parton4.Px();
  m_parton4_py = parton4.Py();
  m_parton4_pz = parton4.Pz();

  m_true_lep1_e = true_lep1.E();
  m_true_lep2_e = true_lep2.E();

  m_true_lep1_px = true_lep1.Px();
  m_true_lep1_py = true_lep1.Py();
  m_true_lep1_pz = true_lep1.Pz();

  m_true_lep2_px = true_lep2.Px();
  m_true_lep2_py = true_lep2.Py();
  m_true_lep2_pz = true_lep2.Pz();

  // -----------   RECO    ---------

  m_jet1_e = jet1.E();
  m_jet2_e = jet2.E();
  m_jet3_e = jet3.E();
  m_jet4_e = jet4.E();

  m_jet1_px = jet1.Px();
  m_jet1_py = jet1.Py();
  m_jet1_pz = jet1.Pz();

  m_jet2_px = jet2.Px();
  m_jet2_py = jet2.Py();
  m_jet2_pz = jet2.Pz();

  m_jet3_px = jet3.Px();
  m_jet3_py = jet3.Py();
  m_jet3_pz = jet3.Pz();

  m_jet4_px = jet4.Px();
  m_jet4_py = jet4.Py();
  m_jet4_pz = jet4.Pz();

  m_lep1_e = lep1.E();
  m_lep2_e = lep2.E();
  
  m_lep1_px = lep1.Px();
  m_lep1_py = lep1.Py();
  m_lep1_pz = lep1.Pz();

  m_lep2_px = lep2.Px();
  m_lep2_py = lep2.Py();
  m_lep2_pz = lep2.Pz();
}

void CompareMEProcessor::processEvent( EVENT::LCEvent *pLCEvent )
{
  this->Clear();
  
  m_nRun = pLCEvent->getRunNumber();
  m_nEvt = pLCEvent->getEventNumber();
  streamlog_out(DEBUG) << "processing event: " << pLCEvent->getEventNumber() << " in run: " << pLCEvent->getRunNumber() << std::endl;

  try {
    // Fetching collections
    LCCollection *inputMCTrueCollection{};
    LCCollection *preselectioncol{};

    streamlog_out(DEBUG) << "  getting true MC collection: " << m_inputMCTrueCollection << std::endl ;
    inputMCTrueCollection = pLCEvent->getCollection( m_inputMCTrueCollection );

    streamlog_out(DEBUG) << "  getting preselection_passed collection: " << m_inputPreSelectionCollection << std::endl ;
    preselectioncol = pLCEvent->getCollection( m_inputPreSelectionCollection );

    // Check for preselection
    m_passed_preselection = preselectioncol->parameters().getIntVal("isPassed");
    if (m_require_presel_pass && !m_passed_preselection)
      return save_evt_with_error_code(ERRORS::PRESELECTION_FAILED_BUT_REQUIRED);

    // Check whether true ZHH or ZZH event
    MCParticle *mcp_H_if_zhh = (MCParticle*) inputMCTrueCollection->getElementAt(10);
    MCParticle *mcp_H_if_zzh = (MCParticle*) inputMCTrueCollection->getElementAt(12);

    // TODO: Check if iterating over inputMCTrueCollection and checking for two H with same parent is "safer"

    m_is_zhh = (mcp_H_if_zhh->getPDG() == 25) && (mcp_H_if_zzh->getPDG() != 25);
    m_is_zzh = (mcp_H_if_zzh->getPDG() == 25) && (mcp_H_if_zhh->getPDG() != 25);

    // Lorentz vectors of final states
    TLorentzVector l1_lortz, l2_lortz;
    TLorentzVector zhh_h1_lortz, zhh_h2_lortz;
    TLorentzVector zzh_z2f1_lortz, zzh_z2f2_lortz, zzh_h_lortz;

    // Save truth info about H/Z-decays
    MCParticle *mcp_l1 = (MCParticle*) inputMCTrueCollection->getElementAt(8); 
    MCParticle *mcp_l2 = (MCParticle*) inputMCTrueCollection->getElementAt(9);

    true_lep1 = v4(mcp_l1);
    true_lep2 = v4(mcp_l2);
    
    MCParticle *mcp_zhh_h1_decay1{}, *mcp_zhh_h1_decay2{}, *mcp_zhh_h2_decay1{}, *mcp_zhh_h2_decay2{};
    MCParticle *mcp_zzh_z2_decay1{}, *mcp_zzh_z2_decay2{}, *mcp_zzh_h1_decay1{}, *mcp_zzh_h1_decay2{};

    // Assign truth-PDGs
    if (m_is_zhh) {
      mcp_zhh_h1_decay1 = (MCParticle*) inputMCTrueCollection->getElementAt(12);
      mcp_zhh_h1_decay2 = (MCParticle*) inputMCTrueCollection->getElementAt(13);
      
      mcp_zhh_h2_decay1 = (MCParticle*) inputMCTrueCollection->getElementAt(14);
      mcp_zhh_h2_decay2 = (MCParticle*) inputMCTrueCollection->getElementAt(15);

      m_true_h1_decay_pdg = abs(mcp_zhh_h1_decay1->getPDG());
      m_true_h2_decay_pdg = abs(mcp_zhh_h2_decay1->getPDG());

      if (m_true_h1_decay_pdg > 0 && m_true_h1_decay_pdg < 7) {
        m_parton1_pdg = m_parton2_pdg = m_true_h1_decay_pdg;

        set_v4(parton1, mcp_zhh_h1_decay1);
        set_v4(parton2, mcp_zhh_h1_decay2);
      }

      if (m_true_h2_decay_pdg > 0 && m_true_h2_decay_pdg < 7) {
        m_parton3_pdg = m_parton4_pdg = m_true_h2_decay_pdg;

        set_v4(parton3, mcp_zhh_h2_decay1);
        set_v4(parton4, mcp_zhh_h2_decay2);
      }
        
    } else if (m_is_zzh) {
      mcp_zzh_z2_decay1 = (MCParticle*) inputMCTrueCollection->getElementAt(10);
      mcp_zzh_z2_decay2 = (MCParticle*) inputMCTrueCollection->getElementAt(11);

      mcp_zzh_h1_decay1 = (MCParticle*) inputMCTrueCollection->getElementAt(13);
      mcp_zzh_h1_decay2 = (MCParticle*) inputMCTrueCollection->getElementAt(14);

      m_true_z2_decay_pdg = abs(mcp_zzh_z2_decay1->getPDG());
      m_true_h1_decay_pdg = abs(mcp_zzh_h1_decay1->getPDG());

      if (m_true_z2_decay_pdg > 0 && m_true_z2_decay_pdg < 7) {
        m_parton1_pdg = m_parton2_pdg = m_true_z2_decay_pdg;

        set_v4(parton1, mcp_zzh_z2_decay1);
        set_v4(parton2, mcp_zzh_z2_decay2);
      }

      if (m_true_h1_decay_pdg > 0 && m_true_h1_decay_pdg < 7) {
        m_parton3_pdg = m_parton4_pdg = m_true_h1_decay_pdg;

        set_v4(parton3, mcp_zzh_h1_decay1);
        set_v4(parton4, mcp_zzh_h1_decay2);
      }
    } else
      return save_evt_with_error_code(ERRORS::UNHANDLED_PROCESS);  
      
    streamlog_out(DEBUG) << "  isZHH " << m_is_zhh << "|" << m_parton1_pdg << ":" << m_parton2_pdg << ":" << m_parton3_pdg << ":" << m_parton4_pdg  << std::endl;

    m_lepton_mode = (m_lepton_mode == -1) ? m_mode : m_lepton_mode;

    // Check ntox pass
    if (m_require_ntox_pass > 0) {
      if (m_require_ntox_pass == 2) {
        if (m_parton1_pdg != m_parton3_pdg || m_parton1_pdg != m_require_ntox) {
          return save_evt_with_error_code(ERRORS::NTOX_FAILED_BUT_REQUIRED);
        }
      }
    }

    // Check TrueJet collection indeed exists (otherwise weird crashs)
    const std::vector<std::string> * colnames = pLCEvent->getCollectionNames();

    // Prepare TrueJet
    TrueJet_Parser * tj = this;
    tj->getall( pLCEvent );

    if (tj->tjcol == NULL)
      return save_evt_with_error_code(ERRORS::TRUEJET_NULL);

    streamlog_out(DEBUG) << " Number of TrueJets is " << tj->njets() << std::endl;
    if ( tj->njets() == 0 )
      return save_evt_with_error_code(ERRORS::TRUEJET_NO_JETS);

    // Save parton-to-true-jet and lepton-to-true-jet-matching, especially for relation between partons and reco jets
    int parton1_jet_i = m_parton1_pdg ? tj->mcpjet(m_is_zhh ? mcp_zhh_h1_decay1 : mcp_zzh_z2_decay1) : -1;
    int parton2_jet_i = m_parton2_pdg ? tj->mcpjet(m_is_zhh ? mcp_zhh_h1_decay2 : mcp_zzh_z2_decay2) : -1;
    int parton3_jet_i = m_parton3_pdg ? tj->mcpjet(m_is_zhh ? mcp_zhh_h2_decay1 : mcp_zzh_h1_decay1) : -1;
    int parton4_jet_i = m_parton4_pdg ? tj->mcpjet(m_is_zhh ? mcp_zhh_h2_decay2 : mcp_zzh_h1_decay2) : -1;

    int lepton1_jet_i = tj->mcpjet(mcp_l1);
    int lepton2_jet_i = tj->mcpjet(mcp_l2);

    streamlog_out(DEBUG) << " TrueJet_IDs|" << parton1_jet_i << ":" << parton2_jet_i << ":" << parton3_jet_i << ":" << parton4_jet_i << ":" << lepton1_jet_i << ":" << lepton2_jet_i << std::endl;

    // A Charged leptons
    if (m_lepton_mode == 0) {
      l1_lortz = v4(mcp_l1);
      l2_lortz = v4(mcp_l2);
    } else if (m_lepton_mode == 1 || m_lepton_mode == 2) {
      ReconstructedParticle *reco_l1{};
      ReconstructedParticle *reco_l2{};

      if (m_lepton_mode == 1) {
        streamlog_out(DEBUG) << "  getting lepton pair: " << m_inputLepPairCollection << std::endl;
        LCCollection *inputLepPair = pLCEvent->getCollection( m_inputLepPairCollection );

        if (inputLepPair->getNumberOfElements() != 2)
          return save_evt_with_error_code(ERRORS::INCOMPLETE_RECO_LEPTON_PAIR);

        reco_l1 = (ReconstructedParticle*) inputLepPair->getElementAt(0);
        reco_l2 = (ReconstructedParticle*) inputLepPair->getElementAt(1);
      } else if (m_lepton_mode == 2) {
        unsigned int n_matching_l_jets = 0;
        unsigned int n_matching_lb_jets = 0;
        for (int i_jet = 0 ; i_jet < tj->njets() ; i_jet++ ) {
          if ( type_jet( i_jet ) == 2 ) { // 2 = leptonic jets
            int jet_pdg = tj->jet(i_jet)->getParticleIDs()[0]->getPDG();
            if ( jet_pdg == m_z1_decay_pdg ) {
              reco_l1 = (ReconstructedParticle*) tj->jet(i_jet);
              n_matching_l_jets++;
            } else if ( jet_pdg == -m_z1_decay_pdg) {
              reco_l2 = (ReconstructedParticle*) tj->jet(i_jet);
              n_matching_lb_jets++;
            }
          }
        }

        if (n_matching_l_jets != n_matching_lb_jets || n_matching_lb_jets != 1)
          return save_evt_with_error_code(ERRORS::INCOMPLETE_TRUEJET_LEPTON_PAIR);
      }

      // Assign final state particles
      l1_lortz = v4(reco_l1);
      l2_lortz = v4(reco_l2);
    } else if (m_lepton_mode == 3) {
      l1_lortz = v4(m_truejet_mode == 0 ? this->p4true(lepton1_jet_i) : this->p4seen(lepton1_jet_i));
      l2_lortz = v4(m_truejet_mode == 0 ? this->p4true(lepton2_jet_i) : this->p4seen(lepton2_jet_i));
    }

    lep1 = l1_lortz;
    lep2 = l2_lortz;

    // B Partons/Jets
    if (m_mode == 0) {
      // True
      MCParticle *mcp_zhh_h1{};
      MCParticle *mcp_zhh_h2{};
      MCParticle *mcp_zzh_h{};

      // Truth mode: Get decay pdg from MCParticle
      m_h1_decay_pdg = m_true_h1_decay_pdg;
      m_h2_decay_pdg = m_true_h2_decay_pdg;
      m_z2_decay_pdg = m_true_z2_decay_pdg;

      if (m_is_zhh) {
        // Assumming ZHH
        mcp_zhh_h1 = (MCParticle*) inputMCTrueCollection->getElementAt(10);
        mcp_zhh_h2 = (MCParticle*) inputMCTrueCollection->getElementAt(11);

        m_z2_decay_mode = getZDecayModeFromPDG(m_h1_decay_pdg);

        zhh_h1_lortz = v4(mcp_zhh_h1);
        zhh_h2_lortz = v4(mcp_zhh_h2);

        m_zhh_is_set = 1;

        // Assuming ZZH (pretend that decay products of H1 are decay products of Z2 in ZZH)          
        zzh_z2f1_lortz = v4(mcp_zhh_h1_decay1);
        zzh_z2f2_lortz = v4(mcp_zhh_h1_decay2);

        if (m_z2_decay_mode > 0) {
          zzh_h_lortz = zhh_h2_lortz; // identify H2 in ZHH as H in ZZH
          m_zzh_is_set = 1;
        }
        
      } else if (m_is_zzh) {
        // Assumming ZZH
        mcp_zzh_h         = (MCParticle*) inputMCTrueCollection->getElementAt(12);

        zzh_z2f1_lortz = v4(mcp_zzh_z2_decay1->getPDG() > 0 ? mcp_zzh_z2_decay1 : mcp_zzh_z2_decay2);
        zzh_z2f2_lortz = v4(mcp_zzh_z2_decay1->getPDG() > 0 ? mcp_zzh_z2_decay2 : mcp_zzh_z2_decay1);

        zzh_h_lortz    = v4(mcp_zzh_h);

        m_z2_decay_mode = getZDecayModeFromPDG(m_z2_decay_pdg);

        m_zzh_is_set = 1;

        // Assuming ZHH
        zhh_h1_lortz = zzh_z2f1_lortz + zzh_z2f2_lortz;
        zhh_h2_lortz = zzh_h_lortz;
        
        m_zhh_is_set = 1;
      }

    } else if (m_mode == 1 || m_mode == 2) {
      // Jet reconstruction mode
      // 1 for RefinedJets
      // 2 for TrueJet


      // For MEM calculation: Use di-jet matching from Misclustering processor
      vector<int> perm_sig; // ZHH
      vector<int> perm_bkg; // ZZH

      vector<ReconstructedParticle*> sig_jets;

      // Check if JetMatching collection exists
      const vector<std::string> *reco_colls = pLCEvent->getCollectionNames();
      if (std::find(reco_colls->begin(), reco_colls->end(), m_inputJetMatchingCollection) == reco_colls->end())
        return save_evt_with_error_code(ERRORS::NO_JET_MATCHING_SIG_COLLECTION);

      // Fetching collections
      streamlog_out(DEBUG) << "  getting JetMatching collection: " << m_inputJetMatchingCollection << std::endl ;

      const EVENT::LCParameters& jm_sig_params = pLCEvent->getCollection( m_inputJetMatchingCollection )->getParameters();

      if (jm_sig_params.getNInt(std::string((m_mode == 1) ? "b2jet2id" : "b2tjet2id")) != 1)
        return save_evt_with_error_code((m_mode == 1) ? ERRORS::INCOMPLETE_RECOJET_COLLECTION : ERRORS::INCOMPLETE_TRUEJET_COLLECTION);

      if (m_mode == 2) {
        // TrueJet mode, using matching by Misclustering processor  

        //DelMe delme(std::bind(&CompareMEProcessor::delall, this));

        perm_sig.push_back(jm_sig_params.getIntVal("b1tjet1id"));
        perm_sig.push_back(jm_sig_params.getIntVal("b1tjet2id"));
        perm_sig.push_back(jm_sig_params.getIntVal("b2tjet1id"));
        perm_sig.push_back(jm_sig_params.getIntVal("b2tjet2id"));

        for (int j = 0; j < 4; j++) {
          sig_jets.push_back((ReconstructedParticle*) tj->jet(perm_sig[j]));
        }

      } else if (m_mode == 1) {
        // From TrueJet to reco matching, and MCParticle to TrueJet matching, we can identify the Reco->MCParticle matching
        // Get TrueJet to reco matching
        std::map<int, int> tj2reco = {
          {jm_sig_params.getIntVal("b1tjet1id"), jm_sig_params.getIntVal("b1jet1id")},
          {jm_sig_params.getIntVal("b1tjet2id"), jm_sig_params.getIntVal("b1jet2id")},
          {jm_sig_params.getIntVal("b2tjet1id"), jm_sig_params.getIntVal("b2jet1id")},
          {jm_sig_params.getIntVal("b2tjet2id"), jm_sig_params.getIntVal("b2jet2id")}
        };

        // Fetching collections
        streamlog_out(DEBUG) << "  getting jet collection: " << m_inputJetCollection << std::endl;
        inputJetCol = pLCEvent->getCollection( m_inputJetCollection );

        // Get jet pairing from Misclustering processor
        perm_sig.push_back(jm_sig_params.getIntVal("b1jet1id"));
        perm_sig.push_back(jm_sig_params.getIntVal("b1jet2id"));
        perm_sig.push_back(jm_sig_params.getIntVal("b2jet1id"));
        perm_sig.push_back(jm_sig_params.getIntVal("b2jet2id"));

        for (int j = 0; j < 4; j++) {
          sig_jets.push_back((ReconstructedParticle*) inputJetCol->getElementAt(perm_sig[j]));
        }

        TString debug_msg;
        debug_msg.Form("TJ->Reco | %i->%i %i | %i->%i %i | %i->%i %i | %i->%i %i",
          jm_sig_params.getIntVal("b1tjet1id"), jm_sig_params.getIntVal("b1jet1id"), int(tj2reco.count(jm_sig_params.getIntVal("b1tjet1id"))),
          jm_sig_params.getIntVal("b1tjet2id"), jm_sig_params.getIntVal("b1jet2id"), int(tj2reco.count(jm_sig_params.getIntVal("b1tjet2id"))),
          jm_sig_params.getIntVal("b2tjet1id"), jm_sig_params.getIntVal("b2jet1id"), int(tj2reco.count(jm_sig_params.getIntVal("b2tjet1id"))),
          jm_sig_params.getIntVal("b2tjet2id"), jm_sig_params.getIntVal("b2jet2id"), int(tj2reco.count(jm_sig_params.getIntVal("b2tjet2id"))) );
          
        streamlog_out(DEBUG) << " TrueJet Matching " << debug_msg << std::endl;

        // Assign et energies if there is a valid mapping in tj2reco
        if (parton1_jet_i > -1 && tj2reco.count(parton1_jet_i)) set_v4(jet1, (ReconstructedParticle*)inputJetCol->getElementAt(tj2reco[parton1_jet_i]));
        if (parton2_jet_i > -1 && tj2reco.count(parton2_jet_i)) set_v4(jet2, (ReconstructedParticle*)inputJetCol->getElementAt(tj2reco[parton2_jet_i]));
        if (parton3_jet_i > -1 && tj2reco.count(parton3_jet_i)) set_v4(jet3, (ReconstructedParticle*)inputJetCol->getElementAt(tj2reco[parton3_jet_i]));
        if (parton4_jet_i > -1 && tj2reco.count(parton4_jet_i)) set_v4(jet4, (ReconstructedParticle*)inputJetCol->getElementAt(tj2reco[parton4_jet_i]));

        debug_msg.Form("Parton:Jet Energies | %.2f:%.2f | %.2f:%.2f | %.2f:%.2f | %.2f:%.2f",
          parton1.E(), jet1.E(),
          parton2.E(), jet2.E(),
          parton3.E(), jet3.E(),
          parton4.E(), jet4.E());
          
        streamlog_out(DEBUG) << debug_msg << std::endl;
      }

      // Save info about regions
      if (jm_sig_params.getNInt("region") == 1)
        m_misclustering_region = jm_sig_params.getIntVal("region");
      if (jm_sig_params.getNInt("region_icns") == 1)
        m_misclustering_region_icns = jm_sig_params.getIntVal("region_icns");

      m_efrac1_reco = (float)jm_sig_params.getDoubleVal("efrac1_reco");
      m_efrac2_reco = (float)jm_sig_params.getDoubleVal("efrac2_reco");
      m_efrac1_true = (float)jm_sig_params.getDoubleVal("efrac1_true");
      m_efrac2_true = (float)jm_sig_params.getDoubleVal("efrac2_true");

      m_efrac1_icn_reco = (float)jm_sig_params.getDoubleVal("efrac1_icn_reco");
      m_efrac2_icn_reco = (float)jm_sig_params.getDoubleVal("efrac2_icn_reco");
      m_efrac1_icn_true = (float)jm_sig_params.getDoubleVal("efrac1_icn_true");
      m_efrac2_icn_true = (float)jm_sig_params.getDoubleVal("efrac2_icn_true");

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
      
      // For now only analyze HH/ZH->bbbarbbbar
      if (!(both_to_b )) // || both_to_c))
        return save_evt_with_error_code(ERRORS::NEITHER_BBBB_NOR_CCCC);
    
      //streamlog_out(DEBUG) << "processEvent : using signal jet pairing " << perm_sig[0] << perm_sig[1] << perm_sig[2] << perm_sig[3] << std::endl;

      // Assuming ZZH
      // b1 = Z
      // b2 = H
      zzh_z2f1_lortz = v4(sig_jets[0]);
      zzh_z2f2_lortz = v4(sig_jets[1]);
      zzh_h_lortz    = v4(sig_jets[2]) + v4(sig_jets[3]);

      m_z2_decay_pdg  = both_to_b ? 5 : 4; // PDGs: bottom->5, charm->4
      m_h1_decay_pdg  = both_to_b ? 5 : 4;
      m_z2_decay_mode = getZDecayModeFromPDG(m_z2_decay_pdg);

      m_zhh_is_set = 1;

      // Assuming ZHH
      // b1,b2 = H
      zhh_h1_lortz = v4(sig_jets[2]) + v4(sig_jets[3]);
      zhh_h2_lortz = v4(sig_jets[0]) + v4(sig_jets[1]);

      m_h2_decay_pdg = both_to_b ? 5 : 4;
      
      m_zzh_is_set = 1;
    } else if (m_mode == 3) {
      // TrueJet only mode, assuming perfect jet-clustering

      // Truth mode: Get decay pdg from MCParticle
      m_h1_decay_pdg = m_true_h1_decay_pdg;
      m_h2_decay_pdg = m_true_h2_decay_pdg;
      m_z2_decay_pdg = m_true_z2_decay_pdg;

      if (!(parton1_jet_i > -1 && parton2_jet_i > -1 && parton3_jet_i > -1 && parton4_jet_i > -1 && lepton1_jet_i > -1 && lepton2_jet_i > -1))
        return save_evt_with_error_code(ERRORS::SOME_TRUEJET_NOT_FOUND);

      // CHEATED: Extract information about decays
      if (m_parton1_pdg != m_parton3_pdg)
        return save_evt_with_error_code(ERRORS::NEITHER_BBBB_NOR_CCCC);

      // Saving jet energies for transfer function (for seen-mode)
      jet1 = v4(m_truejet_mode == 0 ? this->p4true(parton1_jet_i) : this->p4seen(parton1_jet_i));
      jet2 = v4(m_truejet_mode == 0 ? this->p4true(parton2_jet_i) : this->p4seen(parton2_jet_i));
      jet3 = v4(m_truejet_mode == 0 ? this->p4true(parton3_jet_i) : this->p4seen(parton3_jet_i));
      jet4 = v4(m_truejet_mode == 0 ? this->p4true(parton4_jet_i) : this->p4seen(parton4_jet_i));

      // Assuming ZZH
      // b1 = Z
      // b2 = H
      zzh_z2f1_lortz = jet1;
      zzh_z2f2_lortz = jet2;
      zzh_h_lortz    = jet3 + jet4;

      m_z2_decay_pdg  = m_parton1_pdg; // PDGs: bottom->5, charm->4
      m_h1_decay_pdg  = m_parton3_pdg;
      m_z2_decay_mode = getZDecayModeFromPDG(m_z2_decay_pdg);

      m_zhh_is_set = 1;

      // Assuming ZHH
      // b1,b2 = H
      zhh_h1_lortz = jet1 + jet2;
      zhh_h2_lortz = jet3 + jet4;

      m_h2_decay_pdg = m_z2_decay_pdg;
      
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
      TLorentzVector zzh_lortz[5] = {l1_lortz, l2_lortz, zzh_z2f1_lortz, zzh_z2f2_lortz, zzh_h_lortz};
      
      _zzh->SetZDecayMode(m_z1_decay_mode, m_z2_decay_mode);
      _zzh->SetMomentumFinal(zzh_lortz);

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
    }

    transfer_four_momenta();
    m_pTTree->Fill();

    streamlog_out(MESSAGE) << "Processed event | Success | calcZHH:" << m_zhh_is_set << " | calcZZH:" << m_zzh_is_set << std::endl;
    
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

  transfer_four_momenta();
  m_pTTree->Fill();

  streamlog_out(MESSAGE) << "Processed event | Fault | Fault code:" << error_code << std::endl;
}

void CompareMEProcessor::end()
{
  m_pTFile->cd();
  m_pTTree->Write();
  m_pTFile->Close();
  delete m_pTFile;

  //this->delall();
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