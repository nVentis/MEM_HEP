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
 * - Z1DecayMode<int(5)>; defaults to µ+µ- (see ZHH decay modes when running the MEM processor)
 * - HiggsMass<float(125.)>
 * - outputFilename<string("output.root")>
 * - outputTree<string("dataTree")>
*/

#include "CompareTrueMEProcessor.h"
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

CompareTrueMEProcessor aCompareTrueMEProcessor ;

CompareTrueMEProcessor::CompareTrueMEProcessor() :

  Processor("CompareTrueMEProcessor"),
  m_nRun(0),
  m_nEvt(0)

{

	_description = "CompareTrueMEProcessor writes relevant observables to root-file " ;

  registerOutputCollection(LCIO::MCPARTICLE,
				 "MCTrueCollection",
				 "preselection collection",
				 m_inputMCTrueCollection,
				 std::string("preselection")
				 );

	registerProcessorParameter("Z1DecayMode",
        "MEM processor mode of decay of Z1 (=Z for ZHH, Z1 for ZZH)",
        m_Z1DecayMode,
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
        std::string("compare_out.root")
        );

  registerProcessorParameter("outputTree",
        "name of resulting TTree in the output file",
        m_outputTree,
        std::string("dataTree")
        );
}

void CompareTrueMEProcessor::init()
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

  // 1. True
  m_pTTree->Branch("true_is_zhh", &m_true_is_zhh, "true_is_zhh/I");
  m_pTTree->Branch("true_is_zzh", &m_true_is_zzh, "true_is_zzh/I");
  m_pTTree->Branch("true_h1_decay1_pdg", &m_true_h1_decay1_pdg, "true_h1_decay1_pdg/I");

  // 1.a ZHH output
  m_pTTree->Branch("true_zhh_sigma"  , &m_true_zhh_sigma  , "true_zhh_sigma/F");
  m_pTTree->Branch("true_zhh_sigmall", &m_true_zhh_sigmall, "true_zhh_sigmall/F");
  m_pTTree->Branch("true_zhh_sigmalr", &m_true_zhh_sigmalr, "true_zhh_sigmalr/F");
  m_pTTree->Branch("true_zhh_sigmarl", &m_true_zhh_sigmarl, "true_zhh_sigmarl/F");
  m_pTTree->Branch("true_zhh_sigmarr", &m_true_zhh_sigmarr, "true_zhh_sigmarr/F");

  m_pTTree->Branch("true_zhh_mz"  , &m_true_zhh_mz  , "true_zhh_mz/F");
  m_pTTree->Branch("true_zhh_mhh" , &m_true_zhh_mhh , "true_zhh_mhh/F");
  m_pTTree->Branch("true_zhh_mzhh", &m_true_zhh_mzhh, "true_zhh_mzhh/F");

  m_pTTree->Branch("true_zhh_phi" , &m_true_zhh_phi , "true_zhh_phi/F");
  m_pTTree->Branch("true_zhh_phif", &m_true_zhh_phif, "true_zhh_phif/F");
  m_pTTree->Branch("true_zhh_phih", &m_true_zhh_phih, "true_zhh_phih/F");
  
  m_pTTree->Branch("true_zhh_costheta" , &m_true_zhh_costheta , "true_zhh_costheta/F");
  m_pTTree->Branch("true_zhh_costhetaf", &m_true_zhh_costhetaf, "true_zhh_costhetaf/F");
  m_pTTree->Branch("true_zhh_costhetah", &m_true_zhh_costhetah, "true_zhh_costhetah/F");

  // 1.b ZHH input
  m_pTTree->Branch("true_zhh_l1_e" , &m_true_zhh_l1_E , "true_zhh_l1_e/F");
  m_pTTree->Branch("true_zhh_l1_px", &m_true_zhh_l1_px, "true_zhh_l1_px/F");
  m_pTTree->Branch("true_zhh_l1_py", &m_true_zhh_l1_py, "true_zhh_l1_py/F");
  m_pTTree->Branch("true_zhh_l1_pz", &m_true_zhh_l1_pz, "true_zhh_l1_pz/F");

  m_pTTree->Branch("true_zhh_l2_e" , &m_true_zhh_l2_E , "true_zhh_l2_e/F");
  m_pTTree->Branch("true_zhh_l2_px", &m_true_zhh_l2_px, "true_zhh_l2_px/F");
  m_pTTree->Branch("true_zhh_l2_py", &m_true_zhh_l2_py, "true_zhh_l2_py/F");
  m_pTTree->Branch("true_zhh_l2_pz", &m_true_zhh_l2_pz, "true_zhh_l2_pz/F");

  m_pTTree->Branch("true_zhh_h1_e" , &m_true_zhh_h1_E , "true_zhh_h1_e/F");
  m_pTTree->Branch("true_zhh_h1_px", &m_true_zhh_h1_px, "true_zhh_h1_px/F");
  m_pTTree->Branch("true_zhh_h1_py", &m_true_zhh_h1_py, "true_zhh_h1_py/F");
  m_pTTree->Branch("true_zhh_h1_pz", &m_true_zhh_h1_pz, "true_zhh_h1_pz/F");

  m_pTTree->Branch("true_zhh_h2_e" , &m_true_zhh_h2_E , "true_zhh_h2_e/F");
  m_pTTree->Branch("true_zhh_h2_px", &m_true_zhh_h2_px, "true_zhh_h2_px/F");
  m_pTTree->Branch("true_zhh_h2_py", &m_true_zhh_h2_py, "true_zhh_h2_py/F");
  m_pTTree->Branch("true_zhh_h2_pz", &m_true_zhh_h2_pz, "true_zhh_h2_pz/F");


  // 2.a ZZH output
  m_pTTree->Branch("true_zzh_sigma"  , &m_true_zzh_sigma  , "true_zzh_sigma/F");
  
  m_pTTree->Branch("true_zzh_sigmalll", &m_true_zzh_sigmalll, "true_zzh_sigmalll/F");
  m_pTTree->Branch("true_zzh_sigmallr", &m_true_zzh_sigmallr, "true_zzh_sigmallr/F");
  m_pTTree->Branch("true_zzh_sigmalrl", &m_true_zzh_sigmalrl, "true_zzh_sigmalrl/F");
  m_pTTree->Branch("true_zzh_sigmalrr", &m_true_zzh_sigmalrr, "true_zzh_sigmalrr/F");

  m_pTTree->Branch("true_zzh_sigmarrr", &m_true_zzh_sigmalll, "true_zzh_sigmarrr/F");
  m_pTTree->Branch("true_zzh_sigmarrl", &m_true_zzh_sigmarrl, "true_zzh_sigmarrl/F");
  m_pTTree->Branch("true_zzh_sigmarlr", &m_true_zzh_sigmarlr, "true_zzh_sigmarlr/F");
  m_pTTree->Branch("true_zzh_sigmarll", &m_true_zzh_sigmarll, "true_zzh_sigmarll/F");

  m_pTTree->Branch("true_zzh_mz1"  , &m_true_zzh_mz1  , "true_zzh_mz1/F");
  m_pTTree->Branch("true_zzh_mz2"  , &m_true_zzh_mz2  , "true_zzh_mz2/F");
  m_pTTree->Branch("true_zzh_mzz" , &m_true_zzh_mzz , "true_zzh_mzz/F");
  m_pTTree->Branch("true_zzh_mzzh", &m_true_zzh_mzzh, "true_zzh_mzzh/F");

  m_pTTree->Branch("true_zzh_phi" , &m_true_zzh_phi , "true_zzh_phi/F");
  m_pTTree->Branch("true_zzh_phif", &m_true_zzh_phif, "true_zzh_phif/F");
  m_pTTree->Branch("true_zzh_phih", &m_true_zzh_phih, "true_zzh_phih/F");
  
  m_pTTree->Branch("true_zzh_costheta" , &m_true_zzh_costheta , "true_zzh_costheta/F");
  m_pTTree->Branch("true_zzh_costhetaf", &m_true_zzh_costhetaf, "true_zzh_costhetaf/F");
  m_pTTree->Branch("true_zzh_costhetah", &m_true_zzh_costhetah, "true_zzh_costhetah/F");

  // 1.b ZZH input
  m_pTTree->Branch("true_zzh_l1_e" , &m_true_zzh_l1_E , "true_zzh_l1_e/F");
  m_pTTree->Branch("true_zzh_l1_px", &m_true_zzh_l1_px, "true_zzh_l1_px/F");
  m_pTTree->Branch("true_zzh_l1_py", &m_true_zzh_l1_py, "true_zzh_l1_py/F");
  m_pTTree->Branch("true_zzh_l1_pz", &m_true_zzh_l1_pz, "true_zzh_l1_pz/F");

  m_pTTree->Branch("true_zzh_l2_e" , &m_true_zzh_l2_E , "true_zzh_l2_e/F");
  m_pTTree->Branch("true_zzh_l2_px", &m_true_zzh_l2_px, "true_zzh_l2_px/F");
  m_pTTree->Branch("true_zzh_l2_py", &m_true_zzh_l2_py, "true_zzh_l2_py/F");
  m_pTTree->Branch("true_zzh_l2_pz", &m_true_zzh_l2_pz, "true_zzh_l2_pz/F");

  m_pTTree->Branch("true_zzh_h1_e" , &m_true_zzh_h1_E , "true_zzh_h1_e/F");
  m_pTTree->Branch("true_zzh_h1_px", &m_true_zzh_h1_px, "true_zzh_h1_px/F");
  m_pTTree->Branch("true_zzh_h1_py", &m_true_zzh_h1_py, "true_zzh_h1_py/F");
  m_pTTree->Branch("true_zzh_h1_pz", &m_true_zzh_h1_pz, "true_zzh_h1_pz/F");

  m_pTTree->Branch("true_zzh_h2_e" , &m_true_zzh_h2_E , "true_zzh_h2_e/F");
  m_pTTree->Branch("true_zzh_h2_px", &m_true_zzh_h2_px, "true_zzh_h2_px/F");
  m_pTTree->Branch("true_zzh_h2_py", &m_true_zzh_h2_py, "true_zzh_h2_py/F");
  m_pTTree->Branch("true_zzh_h2_pz", &m_true_zzh_h2_pz, "true_zzh_h2_pz/F");

  streamlog_out(DEBUG) << "   init finished  " << std::endl;

  // Print some debug info
  cerr << "\n";
  cerr << "ZHH MEM processor initializing with:\n";
  cerr << "    H_mass: "<< m_Hmass << "\n";
  cerr << "    Z1DecayMode: "<< m_Z1DecayMode << "\n";
  cerr << "    Pol_e: "<< pol_e << "\n";
  cerr << "    Pol_p: "<< pol_p << "\n";
  cerr << "    m_outputTree: " << m_outputTree << "\n";

  _zhh = new LCMEZHH("LCMEZHH", "ZHH", m_Hmass, pol_e, pol_p);
  _zhh->SetZDecayMode(m_Z1DecayMode); // 5 (internal mapping) -> (13) PDG, muon 

  _zzh = new LCMEZZH("LCMEZZH", "ZZH", m_Hmass, pol_e, pol_p);
  // zzh setZDecayMode adjusted every run
}

void CompareTrueMEProcessor::Clear() 
{
  streamlog_out(DEBUG) << "   Clear called  " << std::endl;

  // 1. True
  m_true_is_zhh = 0;
  m_true_is_zzh = 0;
  m_true_h1_decay1_pdg = 0;

  // 2.a ZHH output
  m_true_zhh_sigma     = 0.;
  m_true_zhh_sigmall   = 0.;
  m_true_zhh_sigmalr   = 0.;
  m_true_zhh_sigmarl   = 0.;
  m_true_zhh_sigmarr   = 0.;

  m_true_zhh_mz        = 0.;
  m_true_zhh_mhh       = 0.;
  m_true_zhh_mzhh      = 0.;

  m_true_zhh_phi       = 0.;
  m_true_zhh_phif      = 0.;
  m_true_zhh_phih      = 0.;
  m_true_zhh_costheta  = 0.;
  m_true_zhh_costhetaf = 0.;
  m_true_zhh_costhetah = 0.;

  // 2.b ZHH input
  m_true_zhh_l1_E  = 0.;
  m_true_zhh_l1_px = 0.;
  m_true_zhh_l1_py = 0.;
  m_true_zhh_l1_pz = 0.;

  m_true_zhh_l2_E  = 0.;
  m_true_zhh_l2_px = 0.;
  m_true_zhh_l2_py = 0.;
  m_true_zhh_l2_pz = 0.;

  m_true_zhh_h1_E  = 0.;
  m_true_zhh_h1_px = 0.;
  m_true_zhh_h1_py = 0.;
  m_true_zhh_h1_pz = 0.;


  // 3.a ZZH output
  m_true_zzh_sigma    = 0.;
  m_true_zzh_sigmalll = 0.;
  m_true_zzh_sigmallr = 0.;
  m_true_zzh_sigmalrl = 0.;
  m_true_zzh_sigmalrr = 0.;
  m_true_zzh_sigmarll = 0.;
  m_true_zzh_sigmarlr = 0.;
  m_true_zzh_sigmarrl = 0.;
  m_true_zzh_sigmarrr = 0.;

  m_true_zzh_mz1       = 0.;
  m_true_zzh_mz2       = 0.;
  m_true_zzh_mzz       = 0.;
  m_true_zzh_mzzh      = 0.;
  
  m_true_zzh_phi       = 0.;
  m_true_zzh_phif      = 0.;
  m_true_zzh_phih      = 0.;
  m_true_zzh_costheta  = 0.;
  m_true_zzh_costhetaf = 0.;
  m_true_zzh_costhetah = 0.;

  // 3.b ZZH input
  m_true_zzh_l1_E  = 0.;
  m_true_zzh_l1_px = 0.;
  m_true_zzh_l1_py = 0.;
  m_true_zzh_l1_pz = 0.;

  m_true_zzh_l2_E  = 0.;
  m_true_zzh_l2_px = 0.;
  m_true_zzh_l2_py = 0.;
  m_true_zzh_l2_pz = 0.;

  m_true_zzh_h1_E  = 0.;
  m_true_zzh_h1_px = 0.;
  m_true_zzh_h1_py = 0.;
  m_true_zzh_h1_pz = 0.;

  m_true_zzh_h2_E  = 0.;
  m_true_zzh_h2_px = 0.;
  m_true_zzh_h2_py = 0.;
  m_true_zzh_h2_pz = 0.;
}

void CompareTrueMEProcessor::processRunHeader( LCRunHeader*  /*run*/) { 
  m_nRun++ ;
}

void CompareTrueMEProcessor::processEvent( EVENT::LCEvent *pLCEvent )
{
  this->Clear();
  
  m_nRun = pLCEvent->getRunNumber();
  m_nEvt = pLCEvent->getEventNumber();
  streamlog_out(DEBUG) << "processing event: " << pLCEvent->getEventNumber() << "  in run: " << pLCEvent->getRunNumber() << endl;

  // Helicity combinations
  // ZHH
  Int_t vHelLL[2] = {-1,-1};
  Int_t vHelLR[2] = {-1, 1};
  Int_t vHelRL[2] = { 1,-1};
  Int_t vHelRR[2] = { 1, 1};

  // ZZH
  Int_t vHelLLL[3] = {-1,-1,-1};
  Int_t vHelLLR[3] = {-1,-1,1};
  Int_t vHelLRL[3] = {-1,1,-1};
  Int_t vHelLRR[3] = {-1,1,1};
  Int_t vHelRLL[3] = {1,-1,-1};
  Int_t vHelRLR[3] = {1,-1,1};
  Int_t vHelRRL[3] = {1,1,-1};
  Int_t vHelRRR[3] = {1,1,1};

  LCCollection *inputMCTrueCollection{};

  try {
    // Fetching collections
    streamlog_out(DEBUG) << "        getting true MC collection: " << m_inputMCTrueCollection << std::endl ;
    inputMCTrueCollection = pLCEvent->getCollection( m_inputMCTrueCollection );

    // Check whether true ZHH or ZZH event
    MCParticle *mcPart_H1_if_zhh = dynamic_cast<MCParticle*>(inputMCTrueCollection->getElementAt(10));

    m_true_is_zhh = (mcPart_H1_if_zhh->getPDG() == 25);
    m_true_is_zzh = !m_true_is_zhh;

    // Get particles of final state
    // Same IDs for final Z1 leptons in both true ZZH and ZHH processes
    MCParticle *true_l1 = dynamic_cast<MCParticle*>(inputMCTrueCollection->getElementAt( 8));
    MCParticle *true_l2 = dynamic_cast<MCParticle*>(inputMCTrueCollection->getElementAt( 9));

    TLorentzVector true_l1_lortz = TLorentzVector(true_l1->getMomentum(), true_l1->getEnergy());
    TLorentzVector true_l2_lortz = TLorentzVector(true_l2->getMomentum(), true_l2->getEnergy());

    if (m_true_is_zhh) {
      // Assumming ZHH
      MCParticle *true_zhh_h1 = dynamic_cast<MCParticle*>(inputMCTrueCollection->getElementAt(10));
      MCParticle *true_zhh_h2 = dynamic_cast<MCParticle*>(inputMCTrueCollection->getElementAt(11));

      TLorentzVector true_zhh_h1_lortz = TLorentzVector(true_zhh_h1->getMomentum(), true_zhh_h1->getEnergy());
      TLorentzVector true_zhh_h2_lortz = TLorentzVector(true_zhh_h2->getMomentum(), true_zhh_h2->getEnergy());

      TLorentzVector true_zhh_lortz[4] = {true_l1_lortz, true_l2_lortz, true_zhh_h1_lortz, true_zhh_h2_lortz};

      _zhh->SetMomentumFinal(true_zhh_lortz);

      m_true_zhh_h1_E  = true_zhh_h1->getEnergy();
      m_true_zhh_h1_px = true_zhh_h1->getMomentum()[ 0 ];
      m_true_zhh_h1_py = true_zhh_h1->getMomentum()[ 1 ];
      m_true_zhh_h1_pz = true_zhh_h1->getMomentum()[ 2 ];

      m_true_zhh_h2_E  = true_zhh_h2->getEnergy();
      m_true_zhh_h2_px = true_zhh_h2->getMomentum()[ 0 ];
      m_true_zhh_h2_py = true_zhh_h2->getMomentum()[ 1 ];
      m_true_zhh_h2_pz = true_zhh_h2->getMomentum()[ 2 ];

      // Assuming ZZH
      // Pretend that decay products of H1 are decay products of Z2 in ZZH
      MCParticle *true_zhh_h1_decay1 = dynamic_cast<MCParticle*>(inputMCTrueCollection->getElementAt(12));
      MCParticle *true_zhh_h1_decay2 = dynamic_cast<MCParticle*>(inputMCTrueCollection->getElementAt(13));

      TLorentzVector true_zhh_h1_decay1_lortz = TLorentzVector(true_zhh_h1_decay1->getMomentum(), true_zhh_h1_decay1->getEnergy());
      TLorentzVector true_zhh_h1_decay2_lortz = TLorentzVector(true_zhh_h1_decay2->getMomentum(), true_zhh_h1_decay2->getEnergy());

      m_true_h1_decay1_pdg = true_zhh_h1_decay1->getPDG();

      int z2_decay_mode = getZDecayModeFromPDG(abs(m_true_h1_decay1_pdg));

      if (z2_decay_mode > 0) {
        _zzh->SetZDecayMode(m_Z1DecayMode, z2_decay_mode);

        TLorentzVector true_zzh_lortz[5] = {true_l1_lortz, true_l2_lortz, true_zhh_h1_decay1_lortz, true_zhh_h1_decay2_lortz, true_zhh_h2_lortz};

        _zzh->SetMomentumFinal(true_zzh_lortz);

        // ZZH
        m_true_zzh_mz1   = TMath::Sqrt(_zzh->GetQ2Z1());
        m_true_zzh_mz2   = TMath::Sqrt(_zzh->GetQ2Z2());
        m_true_zzh_mzz  = TMath::Sqrt(_zzh->GetQ2ZZ());
        m_true_zzh_mzzh = TMath::Sqrt(_zzh->GetQ2ZZH());

        m_true_zzh_sigmalll = _zzh->GetMatrixElement2(vHelLLL);
        m_true_zzh_sigmallr = _zzh->GetMatrixElement2(vHelLLR);
        m_true_zzh_sigmalrl = _zzh->GetMatrixElement2(vHelLRL);
        m_true_zzh_sigmalrr = _zzh->GetMatrixElement2(vHelLRR);
        m_true_zzh_sigmarll = _zzh->GetMatrixElement2(vHelRLL);
        m_true_zzh_sigmarlr = _zzh->GetMatrixElement2(vHelRLR);
        m_true_zzh_sigmarrl = _zzh->GetMatrixElement2(vHelRRL);
        m_true_zzh_sigmarrr = _zzh->GetMatrixElement2(vHelRRR);

        m_true_zzh_sigmall = _zzh->GetMatrixElement2(vHelLL);
        m_true_zzh_sigmalr = _zzh->GetMatrixElement2(vHelLR);
        m_true_zzh_sigmarl = _zzh->GetMatrixElement2(vHelRL);
        m_true_zzh_sigmarr = _zzh->GetMatrixElement2(vHelRR);
        m_true_zzh_sigma   = _zzh->GetMatrixElement2();
      }
      
    } else if (m_true_is_zzh) {
      // Get m_Z2DecayMode
      // TODO
    }

    // 1. ZHH
    m_true_zhh_mz   = TMath::Sqrt(_zhh->GetQ2Z());
    m_true_zhh_mhh  = TMath::Sqrt(_zhh->GetQ2HH());
    m_true_zhh_mzhh = TMath::Sqrt(_zhh->GetQ2ZHH());

    m_true_zhh_sigmall = _zhh->GetMatrixElement2(vHelLL);
    m_true_zhh_sigmalr = _zhh->GetMatrixElement2(vHelLR);
    m_true_zhh_sigmarl = _zhh->GetMatrixElement2(vHelRL);
    m_true_zhh_sigmarr = _zhh->GetMatrixElement2(vHelRR);
    m_true_zhh_sigma   = _zhh->GetMatrixElement2();

    m_true_zhh_phi       = _zhh->GetPhi();
    m_true_zhh_phif      = _zhh->GetPhiF();
    m_true_zhh_phih      = _zhh->GetPhiH();
    m_true_zhh_costheta  = _zhh->GetCosTheta();
    m_true_zhh_costhetaf = _zhh->GetCosThetaF();
    m_true_zhh_costhetah = _zhh->GetCosThetaH();

    // 1.b ZHH input
    m_true_zhh_l1_E  = true_l1->getEnergy();
    m_true_zhh_l1_px = true_l1->getMomentum()[ 0 ];
    m_true_zhh_l1_py = true_l1->getMomentum()[ 1 ];
    m_true_zhh_l1_pz = true_l1->getMomentum()[ 2 ];

    m_true_zhh_l2_E  = true_l2->getEnergy();
    m_true_zhh_l2_px = true_l2->getMomentum()[ 0 ];
    m_true_zhh_l2_py = true_l2->getMomentum()[ 1 ];
    m_true_zhh_l2_pz = true_l2->getMomentum()[ 2 ];

    m_pTTree->Fill();
    
  } catch(DataNotAvailableException &e) {
    streamlog_out(MESSAGE) << "processEvent : Input collections not found in event " << m_nEvt << std::endl;
  }

}

void CompareTrueMEProcessor::check()
{
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void CompareTrueMEProcessor::end()
{
  m_pTFile->cd();
  m_pTTree->Write();
  m_pTFile->Close();
  delete m_pTFile;
}

int CompareTrueMEProcessor::getZDecayModeFromPDG(int pdg)
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