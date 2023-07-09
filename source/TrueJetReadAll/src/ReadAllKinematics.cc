#include "ReadAllKinematics.h"
#include <iostream>
#include <vector>
#include <string>
#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>
#include <IMPL/ReconstructedParticleImpl.h>

using namespace lcio ;
using namespace marlin ;
using namespace std ;

template<class T>
TLorentzVector v4(T* p){
  return TLorentzVector( p->getMomentum(), p->getEnergy() );
}

template<class T>
vector<double> fm(T* p){
  return { p->getEnergy(), p->getMomentum()[0], p->getMomentum()[1], p->getMomentum()[2] };
}

template<class T>
double inv_mass(T* p1, T* p2){
  double e = p1->getEnergy()+p2->getEnergy() ;
  double px = p1->getMomentum()[0]+p2->getMomentum()[0];
  double py = p1->getMomentum()[1]+p2->getMomentum()[1];
  double pz = p1->getMomentum()[2]+p2->getMomentum()[2];
  return( sqrt( e*e - px*px - py*py - pz*pz  ) );
}

ReadAllKinematics aReadAllKinematics ;

ReadAllKinematics::ReadAllKinematics() :

  Processor("ReadAllKinematics"),
  n_run(0),
  n_evt(0)
{

	_description = "ReadAllKinematics writes out final state 4 vectors based on MCTruth, TrueJet and reco data" ;

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
				 m_inputJetMatchingSigCollection,
				 std::string("JetMatchingSig")
				 );

  registerInputCollection(LCIO::RECONSTRUCTEDPARTICLE,
				 "JetMatchingBkgCollection",
				 "Name of collection holding parameters about reco/true jet matching for the background hypothesis",
				 m_inputJetMatchingBkgCollection,
				 std::string("JetMatchingBkg")
				 );

  registerProcessorParameter("outputFilename",
        "name of output root file",
        m_outputFile,
        std::string("truejet_all.root")
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

void ReadAllKinematics::init()
{
  streamlog_out(DEBUG) << "   init called  " << std::endl;
  this->Clear();

  // Runtime variables
  n_run = 0;
  n_evt = 0;

  m_pTFile = new TFile(m_outputFile.c_str(), "recreate");
  m_pTTree = new TTree(m_outputTree.c_str(), m_outputTree.c_str());
  m_pTTree->SetDirectory(m_pTFile);

  m_pTTree->Branch("run", &n_run, "run/I");
  m_pTTree->Branch("event", &n_evt, "event/I");
  /**
   * Error code 
   * 
   * 
   * 
  */

  m_pTTree->Branch("passed_preselection", &passed_preselection, "passed_preselection/I");
  m_pTTree->Branch("error_code", &error_code, "error_code/I");

  m_pTTree->Branch("fm_mcp", &fm_mcp);
  m_pTTree->Branch("fm_truejet", &fm_truejet);
  m_pTTree->Branch("fm_recojet", &fm_recojet);

  m_pTTree->Branch("pdgs_mcp", &pdgs_mcp);
  m_pTTree->Branch("pdgs_truejet", &pdgs_truejet);
  m_pTTree->Branch("pdgs_recojet", &pdgs_truejet);

  m_pTTree->Branch("recojet_pid_used", &recojet_pid_used);
  m_pTTree->Branch("truejet_jettypes", &truejet_jettypes);

  streamlog_out(DEBUG) << "   init finished  " << std::endl;
}

void ReadAllKinematics::Clear() 
{
  streamlog_out(DEBUG) << "   Clear called  " << std::endl;

  passed_preselection = 0;
  error_code = 0;

  fm_mcp.clear();
  fm_recojet.clear();
  fm_truejet.clear();

  pdgs_mcp.clear();
  pdgs_recojet.clear();
  pdgs_truejet.clear();

  recojet_pid_used.clear();
  truejet_jettypes.clear();
}

void ReadAllKinematics::processRunHeader( LCRunHeader*  /*run*/) { 
  n_run++ ;
}

void ReadAllKinematics::processEvent( EVENT::LCEvent *pLCEvent )
{
  this->Clear();
  
  n_run = pLCEvent->getRunNumber();
  n_evt = pLCEvent->getEventNumber();

  streamlog_out(DEBUG) << "processing event: " << n_evt << "  in run: " << n_run << endl;

  try {
    // Fetching collections
    LCCollection *inputMCTrueCollection{};
    LCCollection *preselectioncol{};

    streamlog_out(DEBUG) << "        getting true MC collection: " << m_inputMCTrueCollection << std::endl ;
    inputMCTrueCollection = pLCEvent->getCollection( m_inputMCTrueCollection );

    streamlog_out(DEBUG) << "        getting preselection_passed collection: " << m_inputPreSelectionCollection << std::endl ;
    preselectioncol = pLCEvent->getCollection( m_inputPreSelectionCollection );

    streamlog_out(DEBUG) << " getting jet collection: " << m_inputJetCollection << std::endl;
    inputJetCol = pLCEvent->getCollection( m_inputJetCollection );

    // Check for preselection
    passed_preselection = preselectioncol->parameters().getIntVal("isPassed");
    if (require_presel_pass && !passed_preselection)
      return save_evt_with_error_code(ERRORS::PRESELECTION_FAILED_BUT_REQUIRED);

    // Check whether true ZHH or ZZH event
    MCParticle* mcp_H_if_zhh = (MCParticle*) inputMCTrueCollection->getElementAt(10);
    MCParticle* mcp_H_if_zzh = (MCParticle*) inputMCTrueCollection->getElementAt(12);

    bool m_is_zhh = (mcp_H_if_zhh->getPDG() == 25) && (mcp_H_if_zzh->getPDG() != 25);
    bool m_is_zzh = (mcp_H_if_zzh->getPDG() == 25) && (mcp_H_if_zhh->getPDG() != 25);

    // Save MCTruth
    const vector<int> mcp_ids = m_is_zhh ? mcp_ids_zhh : mcp_ids_zzh;

    for (int i = 0; i < mcp_ids.size(); i++) {
      MCParticle* particle = (MCParticle*) inputMCTrueCollection->getElementAt(mcp_ids[i]);
      
      pdgs_mcp.push_back(particle->getPDG());
      fm_mcp.push_back(fm(particle));
    }

    // Save reconstructed jets (if any)
    for (int i = 0; i < inputJetCol->getNumberOfElements(); i++) {
      ReconstructedParticle* particle = (ReconstructedParticle*) inputJetCol->getElementAt(i);

      pdgs_recojet.push_back(particle->getParticleIDUsed()->getPDG());
      recojet_pid_used.push_back(particle->getParticleIDUsed()->getType());
      fm_recojet.push_back(fm(particle));
    }

    // Save TrueJet (if any)
    this->getall( pLCEvent );
    for (int i = 0; this->njets(); i++) {
      ReconstructedParticle* particle = (ReconstructedParticle*) this->jet(i);

      pdgs_truejet.push_back(particle->getParticleIDUsed()->getPDG());
      fm_truejet.push_back(fm(particle));
    }

    m_pTTree->Fill();
    
  } catch(DataNotAvailableException &e) {
    streamlog_out(MESSAGE) << "processEvent : Input collections not found in event " << n_evt << std::endl;
  }

}

void ReadAllKinematics::check()
{
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}

void ReadAllKinematics::save_evt_with_error_code(int error_code)
{
  error_code = error_code;
  m_pTTree->Fill();
}

void ReadAllKinematics::end()
{
  m_pTFile->cd();
  m_pTTree->Write();
  m_pTFile->Close();
  delete m_pTFile;
}
