/**
 * Used to calculate MEs for ZHH and ZZH processes based on pure MC truth data, without TrueJet-jet clustering.
 * TODO: Include flag for using TrueJet output 
 * 
 * @author Bryan Bliewert (TUM/DESY) <bryan.bliewert@tum.de>
 * @date 20.05.2023
 * 
*/

#ifndef CompareTrueMEProcessor_h
#define CompareTrueMEProcessor_h 1

#include "marlin/Processor.h"
#include "IMPL/LCCollectionVec.h"
#include "lcio.h"
#include <string>
#include <TFile.h>
#include <TTree.h>
#include <vector>
#include "TLorentzVector.h"
#include "physsim/LCMEZHH.h"
#include "physsim/LCMEZZH.h"
class TFile;
class TH1F;
class TH1I;
class TH2I;
class TTree;

using namespace lcio ;
using namespace marlin ;
using namespace lcme ;

class CompareTrueMEProcessor : public Processor
{
	public:

		virtual Processor*  newProcessor()
		{
			return new CompareTrueMEProcessor;
		}
		CompareTrueMEProcessor();
		virtual ~CompareTrueMEProcessor() = default;
		CompareTrueMEProcessor(const CompareTrueMEProcessor&) = delete;
		CompareTrueMEProcessor& operator=(const CompareTrueMEProcessor&) = delete;
		virtual void init();
		virtual void Clear();
		virtual void processRunHeader( LCRunHeader*  /*run*/);
		virtual void processEvent( EVENT::LCEvent *pLCEvent );
		virtual void check();
		virtual void end();
		int getZDecayModeFromPDG(int pdg);
		
 protected:
		
		/**
		 * Add the expected output collections
		 */
		
		/** Input collection name.
		 */
		std::string m_inputLepPairCollection{};
		std::string m_inputJetCollection{};
		std::string m_inputHdecayMode{};
		std::string m_inputHiggsPairCollection{};
		std::string m_inputMCTrueCollection{};
		std::string m_outputFile{};
		std::string m_outputTree{};

		int m_nRun;
        int m_nEvt;

		int m_calculation_source{}; // 0 => use MCTruth data; 1 => use reconstructed data (HiggsPair, LeptonPair, HdecayMode, and some jet e.g. RefinedJets)
		float m_Hmass{};

		TFile *m_pTFile{};
        TTree *m_pTTree{};
		lcme::LCMEZHH *_zhh; // ZHH MEM calculator instance
		lcme::LCMEZZH *_zzh; // ZZH MEM calculator instance

		// Control parameters
		int m_zhh_is_set{}; // used internally
		int m_zzh_is_set{}; // used internally

		// True
		int m_z1_decay_mode{}; // as input from parameter; defaults to 5 (mu+mu-)
		int m_z2_decay_mode{}; // when ZZH is assumed; from m_z2_decay1_pdg
		int m_is_zhh{}; // whether true ZHH occured
		int m_is_zzh{}; // whether true ZZH occured
		int m_h1_decay_pdg{}; // if true ZHH occured, abs(PDG) of particle H1 decayed to
		int m_z2_decay_pdg{}; // if true ZZH occured, abs(PDG) of particle Z2 decayed to

		// 1. Assuming ZHH
		// 1.a ZHH output
		float m_zhh_sigma{};
		float m_zhh_sigmall{};
		float m_zhh_sigmalr{};
		float m_zhh_sigmarl{};
		float m_zhh_sigmarr{};

		float m_zhh_mz{};
		float m_zhh_mhh{};
		float m_zhh_mzhh{};

		float m_zhh_phi{};
		float m_zhh_phif{};
		float m_zhh_phih{};

		float m_zhh_costheta{};
		float m_zhh_costhetaf{};
		float m_zhh_costhetah{};

		// 1.b ZHH input
		float m_zhh_l1_E{};
  		float m_zhh_l1_px{};
  		float m_zhh_l1_py{};
  		float m_zhh_l1_pz{};
		
		float m_zhh_l2_E{};
  		float m_zhh_l2_px{};
  		float m_zhh_l2_py{};
  		float m_zhh_l2_pz{};

		float m_zhh_h1_E{};
  		float m_zhh_h1_px{};
  		float m_zhh_h1_py{};
  		float m_zhh_h1_pz{};

		float m_zhh_h2_E{};
  		float m_zhh_h2_px{};
  		float m_zhh_h2_py{};
  		float m_zhh_h2_pz{};

		
		// 1. Assuming ZZH
		// 1.a ZZH output
		float m_zzh_sigma{};
		float m_zzh_sigmalll{};
		float m_zzh_sigmallr{};
		float m_zzh_sigmalrl{};
		float m_zzh_sigmalrr{};
		float m_zzh_sigmarrr{};
		float m_zzh_sigmarrl{};
		float m_zzh_sigmarlr{};
		float m_zzh_sigmarll{};

		float m_zzh_mz1{};
		float m_zzh_mz2{};
		float m_zzh_mzz{};
		float m_zzh_mzzh{};

		float m_zzh_phi{};
		float m_zzh_phiz{};
		float m_zzh_phiz1f{};
		float m_zzh_phiz2f{};

		float m_zzh_costheta{};
		float m_zzh_costhetaz{};
		float m_zzh_costhetaz1f{};
		float m_zzh_costhetaz2f{};

		// 1.b ZZH input
		float m_zzh_l1_E{};
  		float m_zzh_l1_px{};
  		float m_zzh_l1_py{};
  		float m_zzh_l1_pz{};
		
		float m_zzh_l2_E{};
  		float m_zzh_l2_px{};
  		float m_zzh_l2_py{};
  		float m_zzh_l2_pz{};

		float m_zzh_z2f1_E{};
  		float m_zzh_z2f1_px{};
  		float m_zzh_z2f1_py{};
  		float m_zzh_z2f1_pz{};

		float m_zzh_z2f2_E{};
  		float m_zzh_z2f2_px{};
  		float m_zzh_z2f2_py{};
  		float m_zzh_z2f2_pz{};

		float m_zzh_h_E{};
  		float m_zzh_h_px{};
  		float m_zzh_h_py{};
  		float m_zzh_h_pz{};
};

#endif
