/**
 * @author Bryan Bliewert (TUM/DESY) <bryan.bliewert@tum.de>
 * @date 20.05.2023
 * 
 * Based on ZHH .SLCIO files with reconstructed higgs and lepton pairs, calculate the reco and true matrix elements and fetch the corresponding input parameters (all four vectors of truth and reco leptons and higgses).
 * 
 * Parameters:
 * - LepPairCollection<RECONSTRUCTEDPARTICLE>
 * - HiggsPairCollection<RECONSTRUCTEDPARTICLE>
 * - PreSelectionCollection<RECONSTRUCTEDPARTICLE> [why is this a particle and not a simple boolean flag?]
 * - MCTrueCollection<MCPARTICLE>
 * - ZDecayMode<int(5)>; defaults to µ+µ- (see ZHH decay modes when running the MEM processor)
 * - HiggsMass<float(125.)>
 * - outputFilename<string("output.root")>
 * - outputTree<string("dataTree")>
*/

#include "ZHHPostRecoMEProcessor.h"
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

ZHHPostRecoMEProcessor aZHHPostRecoMEProcessor ;

ZHHPostRecoMEProcessor::ZHHPostRecoMEProcessor() :

  Processor("ZHHPostRecoMEProcessor"),
  m_nRun(0),
  m_nEvt(0)

{

	_description = "ZHHPostRecoMEProcessor writes relevant observables to root-file " ;

	registerInputCollection(LCIO::RECONSTRUCTEDPARTICLE,
				"LepPairCollection",
				"Name of reconstructed lepton pair collection",
				m_inputLepPairCollection,
				std::string("LeptonPair")
				);

  registerInputCollection(LCIO::RECONSTRUCTEDPARTICLE,
				"HiggsPairCollection",
				"Name of reconstructed higss pair collection",
				m_inputHiggsPairCollection,
				std::string("HiggsPair")
				);

  registerOutputCollection(LCIO::RECONSTRUCTEDPARTICLE,
				 "PreSelectionCollection",
				 "preselection collection",
				 m_inputPreSelectionCollection,
				 std::string("preselection")
				 );

  registerOutputCollection(LCIO::MCPARTICLE,
				 "MCTrueCollection",
				 "preselection collection",
				 m_inputMCTrueCollection,
				 std::string("preselection")
				 );

	registerProcessorParameter("ZDecayMode",
        "MEM processor mode of ZDecay",
        m_ZDecayMode,
        int(5)
        );

	registerProcessorParameter("HiggsMass",
        "assumed Hmass",
        m_Hmass,
        float(125.)
        );
	
  registerProcessorParameter("outputFilename",
        "name of output root file",
        m_outputFile,
        std::string("output.root")
        );

  registerProcessorParameter("outputTree",
        "name of resulting TTree in the output file",
        m_outputTree,
        std::string("dataTree")
        );
}

void ZHHPostRecoMEProcessor::init()
{
  streamlog_out(DEBUG) << "   init called  " << std::endl;
  this->Clear();

  Double_t pol_e = -1.; // this->parameters()->getFloatVal("beamPol1");
  Double_t pol_p = +1.; // this->parameters()->getFloatVal("beamPol2");

  m_nRun = 0;
  m_nEvt = 0;

  m_pTFile = new TFile(m_outputFile.c_str(),"recreate");
  m_pTTree = new TTree(m_outputTree.c_str(), m_outputTree.c_str());
  m_pTTree->SetDirectory(m_pTFile);

  m_pTTree->Branch("run", &m_nRun, "run/I");
  m_pTTree->Branch("event", &m_nEvt, "event/I");
  m_pTTree->Branch("passed_preselection", &m_passed_preselection, "passed_preselection/I");
  m_pTTree->Branch("passed_validitycheck", &m_passed_validitycheck, "passed_validitycheck/I");
  m_pTTree->Branch("m_true_h1_decay_pdg", &m_true_h1_decay_pdg, "m_true_h1_decay_pdg/I");
  m_pTTree->Branch("m_true_h2_decay_pdg", &m_true_h2_decay_pdg, "m_true_h2_decay_pdg/I");

  // 1. True
  // 1.a ZHH output
  m_pTTree->Branch("true_sigma"  , &m_true_sigma  , "true_sigma/F");
  m_pTTree->Branch("true_sigmall", &m_true_sigmall, "true_sigmall/F");
  m_pTTree->Branch("true_sigmalr", &m_true_sigmalr, "true_sigmalr/F");
  m_pTTree->Branch("true_sigmarl", &m_true_sigmarl, "true_sigmarl/F");
  m_pTTree->Branch("true_sigmarr", &m_true_sigmarr, "true_sigmarr/F");

  m_pTTree->Branch("true_sigma_p1"  , &m_true_sigma_p1  , "true_sigma_p1/F");
  m_pTTree->Branch("true_sigmall_p1", &m_true_sigmall_p1, "true_sigmall_p1/F");
  m_pTTree->Branch("true_sigmalr_p1", &m_true_sigmalr_p1, "true_sigmalr_p1/F");
  m_pTTree->Branch("true_sigmarl_p1", &m_true_sigmarl_p1, "true_sigmarl_p1/F");
  m_pTTree->Branch("true_sigmarr_p1", &m_true_sigmarr_p1, "true_sigmarr_p1/F");

  m_pTTree->Branch("true_sigma_p2"  , &m_true_sigma_p2  , "true_sigma_p2/F");
  m_pTTree->Branch("true_sigmall_p2", &m_true_sigmall_p2, "true_sigmall_p2/F");
  m_pTTree->Branch("true_sigmalr_p2", &m_true_sigmalr_p2, "true_sigmalr_p2/F");
  m_pTTree->Branch("true_sigmarl_p2", &m_true_sigmarl_p2, "true_sigmarl_p2/F");
  m_pTTree->Branch("true_sigmarr_p2", &m_true_sigmarr_p2, "true_sigmarr_p2/F");

  m_pTTree->Branch("true_sigma_p3"  , &m_true_sigma_p3  , "true_sigma_p3/F");
  m_pTTree->Branch("true_sigmall_p3", &m_true_sigmall_p3, "true_sigmall_p3/F");
  m_pTTree->Branch("true_sigmalr_p3", &m_true_sigmalr_p3, "true_sigmalr_p3/F");
  m_pTTree->Branch("true_sigmarl_p3", &m_true_sigmarl_p3, "true_sigmarl_p3/F");
  m_pTTree->Branch("true_sigmarr_p3", &m_true_sigmarr_p3, "true_sigmarr_p3/F");

  m_pTTree->Branch("true_mz"  , &m_true_mz  , "true_mz/F");
  m_pTTree->Branch("true_mhh" , &m_true_mhh , "true_mhh/F");
  m_pTTree->Branch("true_mzhh", &m_true_mzhh, "true_mzhh/F");

  m_pTTree->Branch("true_phi" , &m_true_phi , "true_phi/F");
  m_pTTree->Branch("true_phif", &m_true_phif, "true_phif/F");
  m_pTTree->Branch("true_phih", &m_true_phih, "true_phih/F");
  
  m_pTTree->Branch("true_costheta" , &m_true_costheta , "true_costheta/F");
  m_pTTree->Branch("true_costhetaf", &m_true_costhetaf, "true_costhetaf/F");
  m_pTTree->Branch("true_costhetah", &m_true_costhetah, "true_costhetah/F");

  // 1.b ZHH input
  m_pTTree->Branch("true_l1_e" , &m_true_l1_E , "true_l1_e/F");
  m_pTTree->Branch("true_l1_px", &m_true_l1_px, "true_l1_px/F");
  m_pTTree->Branch("true_l1_py", &m_true_l1_py, "true_l1_py/F");
  m_pTTree->Branch("true_l1_pz", &m_true_l1_pz, "true_l1_pz/F");

  m_pTTree->Branch("true_l2_e" , &m_true_l2_E , "true_l2_e/F");
  m_pTTree->Branch("true_l2_px", &m_true_l2_px, "true_l2_px/F");
  m_pTTree->Branch("true_l2_py", &m_true_l2_py, "true_l2_py/F");
  m_pTTree->Branch("true_l2_pz", &m_true_l2_pz, "true_l2_pz/F");

  m_pTTree->Branch("true_h1_e" , &m_true_h1_E , "true_h1_e/F");
  m_pTTree->Branch("true_h1_px", &m_true_h1_px, "true_h1_px/F");
  m_pTTree->Branch("true_h1_py", &m_true_h1_py, "true_h1_py/F");
  m_pTTree->Branch("true_h1_pz", &m_true_h1_pz, "true_h1_pz/F");

  m_pTTree->Branch("true_h2_e" , &m_true_h2_E , "true_h2_e/F");
  m_pTTree->Branch("true_h2_px", &m_true_h2_px, "true_h2_px/F");
  m_pTTree->Branch("true_h2_py", &m_true_h2_py, "true_h2_py/F");
  m_pTTree->Branch("true_h2_pz", &m_true_h2_pz, "true_h2_pz/F");


  // 2. Reco
  // 2.a ZHH output
  m_pTTree->Branch("reco_sigma"  , &m_reco_sigma  , "reco_sigma/F");
  m_pTTree->Branch("reco_sigmall", &m_reco_sigmall, "reco_sigmall/F");
  m_pTTree->Branch("reco_sigmalr", &m_reco_sigmalr, "reco_sigmalr/F");
  m_pTTree->Branch("reco_sigmarl", &m_reco_sigmarl, "reco_sigmarl/F");
  m_pTTree->Branch("reco_sigmarr", &m_reco_sigmarr, "reco_sigmarr/F");

  m_pTTree->Branch("reco_mz"  , &m_reco_mz  , "reco_mz/F");
  m_pTTree->Branch("reco_mhh" , &m_reco_mhh , "reco_mhh/F");
  m_pTTree->Branch("reco_mzhh", &m_reco_mzhh, "reco_mzhh/F");

  m_pTTree->Branch("reco_phi" , &m_reco_phi , "reco_phi/F");
  m_pTTree->Branch("reco_phif", &m_reco_phif, "reco_phif/F");
  m_pTTree->Branch("reco_phih", &m_reco_phih, "reco_phih/F");

  m_pTTree->Branch("reco_costheta" , &m_reco_costheta , "reco_costheta/F");
  m_pTTree->Branch("reco_costhetaf", &m_reco_costhetaf, "reco_costhetaf/F");
  m_pTTree->Branch("reco_costhetah", &m_reco_costhetah, "reco_costhetah/F");

  // 2.b ZHH input
  m_pTTree->Branch("reco_l1_e" , &m_reco_l1_E , "reco_l1_e/F");
  m_pTTree->Branch("reco_l1_px", &m_reco_l1_px, "reco_l1_px/F");
  m_pTTree->Branch("reco_l1_py", &m_reco_l1_py, "reco_l1_py/F");
  m_pTTree->Branch("reco_l1_pz", &m_reco_l1_pz, "reco_l1_pz/F");

  m_pTTree->Branch("reco_l2_e" , &m_reco_l2_E , "reco_l2_e/F");
  m_pTTree->Branch("reco_l2_px", &m_reco_l2_px, "reco_l2_px/F");
  m_pTTree->Branch("reco_l2_py", &m_reco_l2_py, "reco_l2_py/F");
  m_pTTree->Branch("reco_l2_pz", &m_reco_l2_pz, "reco_l2_pz/F");

  m_pTTree->Branch("reco_h1_e" , &m_reco_h1_E , "reco_h1_e/F");
  m_pTTree->Branch("reco_h1_px", &m_reco_h1_px, "reco_h1_px/F");
  m_pTTree->Branch("reco_h1_py", &m_reco_h1_py, "reco_h1_py/F");
  m_pTTree->Branch("reco_h1_pz", &m_reco_h1_pz, "reco_h1_pz/F");

  m_pTTree->Branch("reco_h2_e" , &m_reco_h2_E , "reco_h2_e/F");
  m_pTTree->Branch("reco_h2_px", &m_reco_h2_px, "reco_h2_px/F");
  m_pTTree->Branch("reco_h2_py", &m_reco_h2_py, "reco_h2_py/F");
  m_pTTree->Branch("reco_h2_pz", &m_reco_h2_pz, "reco_h2_pz/F");

  streamlog_out(DEBUG) << "   init finished  " << std::endl;

  // Print some debug info
  cerr << "\n";
  cerr << "ZHH MEM processor initializing with:\n";
  cerr << "    H_mass: "<< m_Hmass << "\n";
  cerr << "    ZDecayMode: "<< m_ZDecayMode << "\n";
  cerr << "    Pol_e: "<< pol_e << "\n";
  cerr << "    Pol_p: "<< pol_p << "\n";
  cerr << "    m_outputTree: " << m_outputTree << "\n";

  _zhh = new LCMEZHH("LCMEZHH", "ZHH", m_Hmass, pol_e, pol_p);
  _zhh->SetZDecayMode(m_ZDecayMode); // 5 (internal mapping) -> (13) PDG, muon 
}

void ZHHPostRecoMEProcessor::Clear() 
{
  streamlog_out(DEBUG) << "   Clear called  " << std::endl;
  
  m_passed_preselection = 0;
  m_passed_validitycheck = 0;
  m_true_h1_decay_pdg   = 0;
	m_true_h2_decay_pdg   = 0;

  // 1. True
  // 1.a ZHH output
  m_true_sigma_p1     = 0.;
  m_true_sigmall_p1   = 0.;
  m_true_sigmalr_p1   = 0.;
  m_true_sigmarl_p1   = 0.;
  m_true_sigmarr_p1   = 0.;

  m_true_sigma_p2     = 0.;
  m_true_sigmall_p2   = 0.;
  m_true_sigmalr_p2   = 0.;
  m_true_sigmarl_p2   = 0.;
  m_true_sigmarr_p2   = 0.;

  m_true_sigma_p3     = 0.;
  m_true_sigmall_p3   = 0.;
  m_true_sigmalr_p3   = 0.;
  m_true_sigmarl_p3   = 0.;
  m_true_sigmarr_p3   = 0.;

  m_true_sigma     = 0.;
  m_true_sigmall   = 0.;
  m_true_sigmalr   = 0.;
  m_true_sigmarl   = 0.;
  m_true_sigmarr   = 0.;
  m_true_mz        = 0.;
  m_true_mhh       = 0.;
  m_true_mzhh      = 0.;
  m_true_phi       = 0.;
  m_true_phif      = 0.;
  m_true_phih      = 0.;
  m_true_costheta  = 0.;
  m_true_costhetaf = 0.;
  m_true_costhetah = 0.;

  // 2.b ZHH input
  m_true_l1_E  = 0.;
  m_true_l1_px = 0.;
  m_true_l1_py = 0.;
  m_true_l1_pz = 0.;

  m_true_l2_E  = 0.;
  m_true_l2_px = 0.;
  m_true_l2_py = 0.;
  m_true_l2_pz = 0.;

  m_true_h1_E  = 0.;
  m_true_h1_px = 0.;
  m_true_h1_py = 0.;
  m_true_h1_pz = 0.;

  m_true_h2_E  = 0.;
  m_true_h2_px = 0.;
  m_true_h2_py = 0.;
  m_true_h2_pz = 0.;

  
  // 2. Reco
  // 2.a ZHH output
  m_reco_sigma     = 0.;
  m_reco_sigmall   = 0.;
  m_reco_sigmalr   = 0.;
  m_reco_sigmarl   = 0.;
  m_reco_sigmarr   = 0.;
  
  m_reco_mz        = 0.;
  m_reco_mhh       = 0.;
  m_reco_mzhh      = 0.;

  m_reco_phi       = 0.;
  m_reco_phif      = 0.;
  m_reco_phih      = 0.;

  m_reco_costheta  = 0.;
  m_reco_costhetaf = 0.;
  m_reco_costhetah = 0.;

  // 2.b ZHH input
  m_reco_l1_E  = 0.;
  m_reco_l1_px = 0.;
  m_reco_l1_py = 0.;
  m_reco_l1_pz = 0.;

  m_reco_l2_E  = 0.;
  m_reco_l2_px = 0.;
  m_reco_l2_py = 0.;
  m_reco_l2_pz = 0.;

  m_reco_h1_E  = 0.;
  m_reco_h1_px = 0.;
  m_reco_h1_py = 0.;
  m_reco_h1_pz = 0.;

  m_reco_h2_E  = 0.;
  m_reco_h2_px = 0.;
  m_reco_h2_py = 0.;
  m_reco_h2_pz = 0.;
}

void ZHHPostRecoMEProcessor::processRunHeader( LCRunHeader*  /*run*/) { 
  m_nRun++ ;
} 

void ZHHPostRecoMEProcessor::processEvent( EVENT::LCEvent *pLCEvent )
{
  this->Clear();
  
  m_nRun = pLCEvent->getRunNumber();
  m_nEvt = pLCEvent->getEventNumber();
  streamlog_out(DEBUG) << "processing event: " << pLCEvent->getEventNumber() << "  in run: " << pLCEvent->getRunNumber() << endl;

  // Helicity combinations
  Int_t vHelLL[2] = {-1,-1};
  Int_t vHelLR[2] = {-1, 1};
  Int_t vHelRL[2] = { 1,-1};
  Int_t vHelRR[2] = { 1, 1};

  // True
  LCCollection *inputMCTrueCollection{};

  // Reco
  LCCollection *inputLepPairCollection{};
  LCCollection *inputHiggsPairCollection{};
  LCCollection *preselectioncol{};

  try {
    // Fetching collections
    // True
    streamlog_out(DEBUG) << "        getting true MC collection: " << m_inputMCTrueCollection << std::endl ;
    inputMCTrueCollection = pLCEvent->getCollection( m_inputMCTrueCollection );

    // Reco
    streamlog_out(DEBUG) << "        getting lepton pair collection: " << m_inputLepPairCollection << std::endl ;
    inputLepPairCollection = pLCEvent->getCollection( m_inputLepPairCollection );

    streamlog_out(DEBUG) << "        getting higgs pair collection: " << m_inputHiggsPairCollection << std::endl ;
    inputHiggsPairCollection = pLCEvent->getCollection( m_inputHiggsPairCollection );

    streamlog_out(DEBUG) << "        getting preselection_passed collection: " << m_inputPreSelectionCollection << std::endl ;
    preselectioncol = pLCEvent->getCollection( m_inputPreSelectionCollection );

    // Hyperparameters
    // passed_preselection: ~ during (previous) run of PreSelection processor 
    m_passed_preselection = preselectioncol->parameters().getIntVal("isPassed");

    // passed_validitycheck: Exactly two reconstructed leptons and higgses (required for successful reconstruction)
    m_passed_validitycheck = (inputLepPairCollection->getNumberOfElements() == 2 && inputHiggsPairCollection->getNumberOfElements() == 2) ? 1 : 0;

    // Now, fetch ZHH ME input parameters (4 vectors) and calculate ME

    // 1. True
    // 1.a ZHH output
    MCParticle *true_l1 = dynamic_cast<MCParticle*>(inputMCTrueCollection->getElementAt( 8));
    MCParticle *true_l2 = dynamic_cast<MCParticle*>(inputMCTrueCollection->getElementAt( 9));
    
    MCParticle *true_h1 = dynamic_cast<MCParticle*>(inputMCTrueCollection->getElementAt(10));
    MCParticle *true_h2 = dynamic_cast<MCParticle*>(inputMCTrueCollection->getElementAt(11));

    TLorentzVector true_l1_lortz = TLorentzVector(true_l1->getMomentum(), true_l1->getEnergy());
    TLorentzVector true_l2_lortz = TLorentzVector(true_l2->getMomentum(), true_l2->getEnergy());
    TLorentzVector true_h1_lortz = TLorentzVector(true_h1->getMomentum(), true_h1->getEnergy());
    TLorentzVector true_h2_lortz = TLorentzVector(true_h2->getMomentum(), true_h2->getEnergy());

    TLorentzVector true_lortz[4] = {true_l1_lortz, true_l2_lortz, true_h1_lortz, true_h2_lortz};
    TLorentzVector true_lortz_p1[4] = {true_l1_lortz, true_l2_lortz, true_h2_lortz, true_h1_lortz}; // H1 <-> H2
    TLorentzVector true_lortz_p2[4] = {true_l2_lortz, true_l1_lortz, true_h1_lortz, true_h2_lortz}; // L1 <-> L2
    TLorentzVector true_lortz_p3[4] = {true_l2_lortz, true_l1_lortz, true_h2_lortz, true_h1_lortz}; // L1 <-> L2 AND H1 <-> H2

    _zhh->SetMomentumFinal(true_lortz_p1);

    m_true_sigmall_p1 = _zhh->GetMatrixElement2(vHelLL);
    m_true_sigmalr_p1 = _zhh->GetMatrixElement2(vHelLR);
    m_true_sigmarl_p1 = _zhh->GetMatrixElement2(vHelRL);
    m_true_sigmarr_p1 = _zhh->GetMatrixElement2(vHelRR);
    m_true_sigma_p1   = _zhh->GetMatrixElement2();

    _zhh->SetMomentumFinal(true_lortz_p2);

    m_true_sigmall_p2 = _zhh->GetMatrixElement2(vHelLL);
    m_true_sigmalr_p2 = _zhh->GetMatrixElement2(vHelLR);
    m_true_sigmarl_p2 = _zhh->GetMatrixElement2(vHelRL);
    m_true_sigmarr_p2 = _zhh->GetMatrixElement2(vHelRR);
    m_true_sigma_p2   = _zhh->GetMatrixElement2();

    _zhh->SetMomentumFinal(true_lortz_p3);

    m_true_sigmall_p3 = _zhh->GetMatrixElement2(vHelLL);
    m_true_sigmalr_p3 = _zhh->GetMatrixElement2(vHelLR);
    m_true_sigmarl_p3 = _zhh->GetMatrixElement2(vHelRL);
    m_true_sigmarr_p3 = _zhh->GetMatrixElement2(vHelRR);
    m_true_sigma_p3   = _zhh->GetMatrixElement2();
    
    _zhh->SetMomentumFinal(true_lortz);

    m_true_sigmall = _zhh->GetMatrixElement2(vHelLL);
    m_true_sigmalr = _zhh->GetMatrixElement2(vHelLR);
    m_true_sigmarl = _zhh->GetMatrixElement2(vHelRL);
    m_true_sigmarr = _zhh->GetMatrixElement2(vHelRR);
    m_true_sigma   = _zhh->GetMatrixElement2();

    m_true_mz   = TMath::Sqrt(_zhh->GetQ2Z());
    m_true_mhh  = TMath::Sqrt(_zhh->GetQ2HH());
    m_true_mzhh = TMath::Sqrt(_zhh->GetQ2ZHH());

    m_true_phi       = _zhh->GetPhi();
    m_true_phif      = _zhh->GetPhiF();
    m_true_phih      = _zhh->GetPhiH();
    m_true_costheta  = _zhh->GetCosTheta();
    m_true_costhetaf = _zhh->GetCosThetaF();
    m_true_costhetah = _zhh->GetCosThetaH();

    MCParticle *mcPart_h1_decay = dynamic_cast<MCParticle*>(inputMCTrueCollection->getElementAt(12));
    m_true_h1_decay_pdg = abs(mcPart_h1_decay->getPDG());

    MCParticle *mcPart_h2_decay = dynamic_cast<MCParticle*>(inputMCTrueCollection->getElementAt(14));
    m_true_h2_decay_pdg = abs(mcPart_h2_decay->getPDG());

    // 1.b ZHH input
    m_true_l1_E  = true_l1->getEnergy();
    m_true_l1_px = true_l1->getMomentum()[ 0 ];
    m_true_l1_py = true_l1->getMomentum()[ 1 ];
    m_true_l1_pz = true_l1->getMomentum()[ 2 ];

    m_true_l2_E  = true_l2->getEnergy();
    m_true_l2_px = true_l2->getMomentum()[ 0 ];
    m_true_l2_py = true_l2->getMomentum()[ 1 ];
    m_true_l2_pz = true_l2->getMomentum()[ 2 ];

    m_true_h1_E  = true_h1->getEnergy();
    m_true_h1_px = true_h1->getMomentum()[ 0 ];
    m_true_h1_py = true_h1->getMomentum()[ 1 ];
    m_true_h1_pz = true_h1->getMomentum()[ 2 ];

    m_true_h2_E  = true_h2->getEnergy();
    m_true_h2_px = true_h2->getMomentum()[ 0 ];
    m_true_h2_py = true_h2->getMomentum()[ 1 ];
    m_true_h2_pz = true_h2->getMomentum()[ 2 ];


    // 2. Reco
    if (m_passed_validitycheck) {
      ReconstructedParticle* reco_l1 = static_cast<ReconstructedParticle*>(inputLepPairCollection->getElementAt(0));
      ReconstructedParticle* reco_l2 = static_cast<ReconstructedParticle*>(inputLepPairCollection->getElementAt(1));
      
      ReconstructedParticle* reco_h1 = static_cast<ReconstructedParticle*>(inputHiggsPairCollection->getElementAt(0));
      ReconstructedParticle* reco_h2 = static_cast<ReconstructedParticle*>(inputHiggsPairCollection->getElementAt(1));

      //TLorentzVector l1_reco_lortz = TLorentzVector( l1_reco->getMomentum()[ 0 ] , l1_reco->getMomentum()[ 1 ] , l1_reco->getMomentum()[ 2 ] , l1_reco->getEnergy() );
      //TLorentzVector l2_reco_lortz = TLorentzVector( l2_reco->getMomentum()[ 0 ] , l2_reco->getMomentum()[ 1 ] , l2_reco->getMomentum()[ 2 ] , l2_reco->getEnergy() );
      
      //TLorentzVector h1_reco_lortz = TLorentzVector( h1_reco->getMomentum()[ 0 ] , h1_reco->getMomentum()[ 1 ] , h1_reco->getMomentum()[ 2 ] , h1_reco->getEnergy() );
      //TLorentzVector h2_reco_lortz = TLorentzVector( h2_reco->getMomentum()[ 0 ] , h2_reco->getMomentum()[ 1 ] , h2_reco->getMomentum()[ 2 ] , h2_reco->getEnergy() );

      TLorentzVector reco_l1_lortz = TLorentzVector(reco_l1->getMomentum(), reco_l1->getEnergy());
      TLorentzVector reco_l2_lortz = TLorentzVector(reco_l2->getMomentum(), reco_l2->getEnergy());
      TLorentzVector reco_h1_lortz = TLorentzVector(reco_h1->getMomentum(), reco_h1->getEnergy());
      TLorentzVector reco_h2_lortz = TLorentzVector(reco_h2->getMomentum(), reco_h2->getEnergy());

      TLorentzVector reco_lortz[4] = {reco_l1_lortz, reco_l2_lortz, reco_h1_lortz, reco_h2_lortz};

      _zhh->SetMomentumFinal(reco_lortz);

      // 2.a ZHH output
      m_reco_sigmall   = _zhh->GetMatrixElement2(vHelLL);
      m_reco_sigmalr   = _zhh->GetMatrixElement2(vHelLR);
      m_reco_sigmarl   = _zhh->GetMatrixElement2(vHelRL);
      m_reco_sigmarr   = _zhh->GetMatrixElement2(vHelRR);
      m_reco_sigma     = _zhh->GetMatrixElement2();

      m_reco_mz        = TMath::Sqrt(_zhh->GetQ2Z());
      m_reco_mhh       = TMath::Sqrt(_zhh->GetQ2HH());
      m_reco_mzhh      = TMath::Sqrt(_zhh->GetQ2ZHH());

      m_reco_phi       = _zhh->GetPhi();
      m_reco_phif      = _zhh->GetPhiF();
      m_reco_phih      = _zhh->GetPhiH();
      m_reco_costheta  = _zhh->GetCosTheta();
      m_reco_costhetaf = _zhh->GetCosThetaF();
      m_reco_costhetah = _zhh->GetCosThetaH();

      // 2.b ZHH input
      m_reco_l1_E  = reco_l1->getEnergy();
      m_reco_l1_px = reco_l1->getMomentum()[ 0 ];
      m_reco_l1_py = reco_l1->getMomentum()[ 1 ];
      m_reco_l1_pz = reco_l1->getMomentum()[ 2 ];

      m_reco_l2_E  = reco_l2->getEnergy();
      m_reco_l2_px = reco_l2->getMomentum()[ 0 ];
      m_reco_l2_py = reco_l2->getMomentum()[ 1 ];
      m_reco_l2_pz = reco_l2->getMomentum()[ 2 ];

      m_reco_h1_E  = reco_h1->getEnergy();
      m_reco_h1_px = reco_h1->getMomentum()[ 0 ];
      m_reco_h1_py = reco_h1->getMomentum()[ 1 ];
      m_reco_h1_pz = reco_h1->getMomentum()[ 2 ];

      m_reco_h2_E  = reco_h2->getEnergy();
      m_reco_h2_px = reco_h2->getMomentum()[ 0 ];
      m_reco_h2_py = reco_h2->getMomentum()[ 1 ];
      m_reco_h2_pz = reco_h2->getMomentum()[ 2 ];
    }

    m_pTTree->Fill();
    
  } catch(DataNotAvailableException &e) {
    streamlog_out(MESSAGE) << "processEvent : Input collections not found in event " << m_nEvt << std::endl;
  }

}

void ZHHPostRecoMEProcessor::check()
{
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void ZHHPostRecoMEProcessor::end()
{
  m_pTFile->cd();
  m_pTTree->Write();
  m_pTFile->Close();
  delete m_pTFile;
}
