#include "SaveFinalKinematics.h"
#include <iostream>
#include <vector>
#include <string>
#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>
#include <IMPL/ReconstructedParticleImpl.h>
#include <UTIL/PIDHandler.h>
#include <TMath.h>

using namespace lcio ;
using namespace marlin ;

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

SaveFinalKinematics aSaveFinalKinematics ;

SaveFinalKinematics::SaveFinalKinematics() :

  Processor("SaveFinalKinematics"),
  n_run(0),
  n_evt(0)
{

	_description = "SaveFinalKinematics writes out final state 4 vectors based on MCTruth, TrueJet and reco data" ;

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
				 "RecoParticleCollection",
				 "collection of particle flow objects",
				 m_inputRecoParticleCollection,
				 std::string("PandoraPFOs")
				 );

  registerProcessorParameter("recoPfoAlgoType",
        "which reconstruction algo to use for pfo pid",
        m_recoPfoAlgoType,
        std::string("LikelihoodPID") // as of now, only lcfiplus implemented
        );

  registerProcessorParameter("recoJetAlgoType",
        "which reconstruction type to use for jet pid",
        m_recoJetAlgoType,
        std::string("lcfiplus") // as of now, only lcfiplus implemented
        );

  registerProcessorParameter("minBLikeliness",
        "threshold necessary for BTag value",
        m_minBLikeliness,
        float(0.2)
        );

  registerProcessorParameter("minCLikeliness",
        "threshold necessary for CTag value",
        m_minCLikeliness,
        float(0.2)
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

void SaveFinalKinematics::init()
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
  m_pTTree->Branch("fm_recopfo", &fm_recopfo);
  m_pTTree->Branch("fm_recojet", &fm_recojet);
  m_pTTree->Branch("fm_truejet", &fm_truejet);

  m_pTTree->Branch("pdg_mcp", &pdg_mcp);
  m_pTTree->Branch("pdg_recopfo", &pdg_recopfo);
  m_pTTree->Branch("pdg_recojet", &pdg_recojet);
  m_pTTree->Branch("pdg_truejet", &pdg_truejet);

  m_pTTree->Branch("charge_mcp", &charge_mcp);
  m_pTTree->Branch("charge_recopfo", &charge_recopfo);
  m_pTTree->Branch("charge_recojet", &charge_recojet);
  m_pTTree->Branch("charge_truejet", &charge_truejet);


  m_pTTree->Branch("truejet_jettypes", &truejet_jettypes);

  debug_print = 1;

  streamlog_out(DEBUG) << "   init finished  " << std::endl;
}

void SaveFinalKinematics::Clear() 
{
  streamlog_out(DEBUG) << "   Clear called  " << std::endl;

  passed_preselection = 0;
  error_code = 0;

  fm_mcp.clear();
  fm_recojet.clear();
  fm_recopfo.clear();
  fm_truejet.clear();

  pdg_mcp.clear();
  pdg_recojet.clear();
  pdg_recopfo.clear();
  pdg_truejet.clear();

  charge_mcp.clear();
	charge_recojet.clear();
	charge_recopfo.clear();
	charge_truejet.clear();

  truejet_jettypes.clear();
}

void SaveFinalKinematics::processRunHeader( LCRunHeader*  /*run*/) { 
  n_run++ ;
}

void SaveFinalKinematics::processEvent( EVENT::LCEvent *pLCEvent )
{
  this->Clear();
  
  n_run = pLCEvent->getRunNumber();
  n_evt = pLCEvent->getEventNumber();

  streamlog_out(DEBUG) << "processing event: " << n_evt << "  in run: " << n_run << std::endl;

  try {
    // Check for preselection
    if (require_presel_pass) {
      streamlog_out(DEBUG) << "        getting preselection_passed collection: " << m_inputPreSelectionCollection << std::endl ;
      LCCollection *preselectioncol = pLCEvent->getCollection( m_inputPreSelectionCollection );

      if (!preselectioncol->parameters().getIntVal("isPassed"))
        return save_evt_with_error_code(ERRORS::PRESELECTION_FAILED_BUT_REQUIRED);
    }

    // Fetching collections
    LCCollection *inputMCTrueCollection{};
    LCCollection *inputPfoCol{};
    LCCollection *inputJetCol{};

    streamlog_out(DEBUG) << "        getting true MC collection: " << m_inputMCTrueCollection << std::endl ;
    inputMCTrueCollection = pLCEvent->getCollection( m_inputMCTrueCollection );

    streamlog_out(DEBUG) << "        getting particle flow object collection: " << m_inputRecoParticleCollection << std::endl;
    inputPfoCol = pLCEvent->getCollection( m_inputRecoParticleCollection );

    streamlog_out(DEBUG) << "        getting jet collection: " << m_inputJetCollection << std::endl;
    inputJetCol = pLCEvent->getCollection( m_inputJetCollection );

    // Save MCTruth
    MCParticle* mcp_H_if_zhh = (MCParticle*) inputMCTrueCollection->getElementAt(10);
    MCParticle* mcp_H_if_zzh = (MCParticle*) inputMCTrueCollection->getElementAt(12);

    bool m_is_zhh = (mcp_H_if_zhh->getPDG() == 25) && (mcp_H_if_zzh->getPDG() != 25);
    bool m_is_zzh = (mcp_H_if_zzh->getPDG() == 25) && (mcp_H_if_zhh->getPDG() != 25);

    // Check whether true ZHH or ZZH event
    if (m_is_zhh || m_is_zzh) {
      const vector<int> mcp_ids = m_is_zhh ? mcp_ids_zhh : mcp_ids_zzh;

      for (int i = 0; i < (int)mcp_ids.size(); i++) {
        MCParticle* particle = (MCParticle*) inputMCTrueCollection->getElementAt(mcp_ids[i]);

        pdg_mcp.push_back(particle->getPDG());
        fm_mcp.push_back(fm(particle));
        charge_mcp.push_back(particle->getCharge());
      }
    }

    // Save PandoraPFOs
    PIDHandler pidh_pfo ( inputPfoCol );
    int algo_pfo = pidh_pfo.getAlgorithmID(m_recoPfoAlgoType);
    for (int i = 0; i < inputPfoCol->getNumberOfElements(); i++) {
      ReconstructedParticle* particle = (ReconstructedParticle*) inputPfoCol->getElementAt(i);
      float charge = particle->getCharge();

      const ParticleID &pid = pidh_pfo.getParticleID(particle, algo_pfo);
      int pdg = TMath::Sign(pid.getPDG(), charge);

      pdg_recopfo.push_back(pdg);

      std::cerr << "Particle ["  << i << "] PDG=[" << pdg << "]" << std::endl;

      fm_recopfo.push_back(fm(particle));
      charge_recopfo.push_back(charge);
    }

    // Save reconstructed jets (if any)
    // As of now, with only using lcfiplus, only BTagging and CTagging are supported
    PIDHandler pidh_jet ( inputJetCol );
    for (int i = 0; i < inputJetCol->getNumberOfElements(); i++) {
      ReconstructedParticle* particle = (ReconstructedParticle*) inputJetCol->getElementAt(i);
      float charge = particle->getCharge();

      //ParticleIDImpl* inPID = dynamic_cast<ParticleIDImpl*>(outputPFO->getParticleIDs()[j]);

      if (debug_print) {
        const EVENT::LCParameters& jet_params = pLCEvent->getCollection("RefinedJets")->getParameters();

        ParticleIDVec pids = particle->getParticleIDs();

        vector<std::string> pid_algo_names;
        jet_params.getStringVals("PIDAlgorithmTypeName", pid_algo_names);

        for (int j=0; j < (int)pids.size(); j++) {
          std::string pid_algo_name = pid_algo_names[j];

          vector<std::string> pid_param_names;
          jet_params.getStringVals("ParameterNames_"+pid_algo_name, pid_param_names);

          ParticleID *pid = pids[j];
          FloatVec params = pid->getParameters();

          std::cerr << "PIDAlgo ["  << j <<":" << pid_algo_name <<"]" << std::endl;
          for (int k=0; k < (int)pid_param_names.size(); k++) {
            std::cerr << "Parameter ["  << pid_param_names[k] <<"]=[" << params[k] << "]" << std::endl;
          }
          std::cerr << std::endl;
        }
      }

      if (m_recoJetAlgoType == "lcfiplus") {
        int algo_jet = pidh_jet.getAlgorithmID("lcfiplus");
        int ibtag = pidh_jet.getParameterIndex(algo_jet, "BTag");
        int ictag = pidh_jet.getParameterIndex(algo_jet, "CTag");

        const ParticleID &pid = pidh_jet.getParticleID(particle, algo_jet);

        float bLikeliness = pid.getParameters()[ibtag];
        float cLikeliness = pid.getParameters()[ictag];

        //std::cerr << "RecoPID of " << j << "(" << pids.size() << "):" << pid->getAlgorithmType() << "=" << pid->getPDG() << std::endl;
        std::cerr << "Particle ["  << i << "] BTag=[" << bLikeliness << "] CTag=[" << cLikeliness << "]" << std::endl;

        if (bLikeliness >= m_minBLikeliness && bLikeliness > cLikeliness)
          pdg_recojet.push_back(TMath::Sign(5, charge));
        else if (cLikeliness >= m_minCLikeliness && cLikeliness > bLikeliness)
          pdg_recojet.push_back(TMath::Sign(4, charge));
        else
          pdg_recojet.push_back(0);

      } else
        throw UnknownAlgorithm(m_recoJetAlgoType);

      fm_recojet.push_back(fm(particle));
      charge_recojet.push_back(particle->getCharge());
    }

    if (debug_print) {
      debug_print=0;
    }

    // Save TrueJet (if any)
    this->getall( pLCEvent );
    for (int i = 0; i < this->njets(); i++) {
      ReconstructedParticle* particle = (ReconstructedParticle*) this->jet(i);
      ParticleID* pid = particle->getParticleIDs()[0];

      pdg_truejet.push_back(pid->getPDG());
      fm_truejet.push_back(fm(particle));
      charge_truejet.push_back(particle->getCharge());

      truejet_jettypes.push_back(type_jet(i));
    }

    m_pTTree->Fill();
    
  } catch(DataNotAvailableException &e) {
    streamlog_out(MESSAGE) << "processEvent : Input collections not found in event " << n_evt << std::endl;
  }

}

void SaveFinalKinematics::check()
{
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}

void SaveFinalKinematics::save_evt_with_error_code(int i_error_code)
{
  error_code = i_error_code;
  m_pTTree->Fill();
}

void SaveFinalKinematics::end()
{
  m_pTFile->cd();
  m_pTTree->Write();
  m_pTFile->Close();
  delete m_pTFile;
}
