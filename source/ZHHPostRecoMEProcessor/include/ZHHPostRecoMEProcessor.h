/**
 * @author Bryan Bliewert (TUM/DESY) <bryan.bliewert@tum.de>
 * @date 20.05.2023
 * 
*/

#ifndef ZHHPostRecoMEProcessor_h
#define ZHHPostRecoMEProcessor_h 1

#include "marlin/Processor.h"
#include "IMPL/LCCollectionVec.h"
#include "lcio.h"
#include <string>
#include <TFile.h>
#include <TTree.h>
#include <vector>
#include "TLorentzVector.h"
#include "physsim/LCMEZHH.h"
class TFile;
class TH1F;
class TH1I;
class TH2I;
class TTree;

using namespace lcio ;
using namespace marlin ;
using namespace lcme ;

class ZHHPostRecoMEProcessor : public Processor
{
	public:

		virtual Processor*  newProcessor()
		{
			return new ZHHPostRecoMEProcessor;
		}
		ZHHPostRecoMEProcessor();
		virtual ~ZHHPostRecoMEProcessor() = default;
		ZHHPostRecoMEProcessor(const ZHHPostRecoMEProcessor&) = delete;
		ZHHPostRecoMEProcessor& operator=(const ZHHPostRecoMEProcessor&) = delete;
		virtual void init();
		virtual void Clear();
		virtual void processRunHeader( LCRunHeader*  /*run*/);
		virtual void processEvent( EVENT::LCEvent *pLCEvent );
		virtual void check();
		virtual void end();
		
 protected:
		
		/**
		 * Add the expected output collections
		 */
		
		/** Input collection name.
		 */
		std::string m_inputLepPairCollection{};
		std::string m_inputHiggsPairCollection{};
		std::string m_inputPreSelectionCollection{};
		std::string m_inputMCTrueCollection{};
		std::string m_outputFile{};
		std::string m_outputTree{};

		int m_nRun;
        int m_nEvt;
		int m_ZDecayMode{};
		float m_Hmass{};

		TFile *m_pTFile{};
        TTree *m_pTTree{};
		lcme::LCMEZHH *_zhh; // ZHH MEM calculator instance

		// Refreshed every run
		int m_passed_preselection{};
		int m_passed_validitycheck{};
		int m_true_h1_decay_pdg{};
		int m_true_h2_decay_pdg{};

		// 1. True
		// 1.a ZHH output
		float m_true_sigma{};
		float m_true_sigmall{};
		float m_true_sigmalr{};
		float m_true_sigmarl{};
		float m_true_sigmarr{};

		float m_true_sigma_p1{};
		float m_true_sigmall_p1{};
		float m_true_sigmalr_p1{};
		float m_true_sigmarl_p1{};
		float m_true_sigmarr_p1{};

		float m_true_sigma_p2{};
		float m_true_sigmall_p2{};
		float m_true_sigmalr_p2{};
		float m_true_sigmarl_p2{};
		float m_true_sigmarr_p2{};

		float m_true_sigma_p3{};
		float m_true_sigmall_p3{};
		float m_true_sigmalr_p3{};
		float m_true_sigmarl_p3{};
		float m_true_sigmarr_p3{};

		float m_true_mz{};
		float m_true_mhh{};
		float m_true_mzhh{};

		float m_true_phi{};
		float m_true_phif{};
		float m_true_phih{};

		float m_true_costheta{};
		float m_true_costhetaf{};
		float m_true_costhetah{};

		// 1.b ZHH input
		float m_true_l1_E{};
  		float m_true_l1_px{};
  		float m_true_l1_py{};
  		float m_true_l1_pz{};
		
		float m_true_l2_E{};
  		float m_true_l2_px{};
  		float m_true_l2_py{};
  		float m_true_l2_pz{};

		float m_true_h1_E{};
  		float m_true_h1_px{};
  		float m_true_h1_py{};
  		float m_true_h1_pz{};

		float m_true_h2_E{};
  		float m_true_h2_px{};
  		float m_true_h2_py{};
  		float m_true_h2_pz{};


		// 2. Reco
		// 2.a ZHH output
		float m_reco_sigma{};
		float m_reco_sigmall{};
		float m_reco_sigmalr{};
		float m_reco_sigmarl{};
		float m_reco_sigmarr{};

		float m_reco_mz{};
		float m_reco_mhh{};
		float m_reco_mzhh{};

		float m_reco_phi{};
		float m_reco_phif{};
		float m_reco_phih{};

		float m_reco_costheta{};
		float m_reco_costhetaf{};
		float m_reco_costhetah{};

		// 2.b ZHH input
		float m_reco_l1_E{};
  		float m_reco_l1_px{};
  		float m_reco_l1_py{};
  		float m_reco_l1_pz{};
		
		float m_reco_l2_E{};
  		float m_reco_l2_px{};
  		float m_reco_l2_py{};
  		float m_reco_l2_pz{};

		float m_reco_h1_E{};
  		float m_reco_h1_px{};
  		float m_reco_h1_py{};
  		float m_reco_h1_pz{};

		float m_reco_h2_E{};
  		float m_reco_h2_px{};
  		float m_reco_h2_py{};
  		float m_reco_h2_pz{};
};

#endif
