#ifndef SaveFinalKinematics_h
#define SaveFinalKinematics_h 1

#include "marlin/Processor.h"
#include "IMPL/LCCollectionVec.h"
#include "lcio.h"
#include <string>
#include <TFile.h>
#include <TTree.h>
#include <vector>
#include "TLorentzVector.h"
#include "TrueJet_Parser.h"
class TFile;
class TH1F;
class TH1I;
class TH2I;
class TTree;

using namespace lcio;
using namespace marlin;
using namespace std;

enum ERRORS : unsigned int {
  PRESELECTION_FAILED_BUT_REQUIRED = 1,
  INCOMPLETE_RECO_LEPTON_PAIR = 101,
  INCOMPLETE_RECOJET_COLLECTION = 102,

  INCOMPLETE_TRUEJET_LEPTON_PAIR = 201,
  INCOMPLETE_TRUEJET_COLLECTION = 202,
  
  NEITHER_BBBB_NOR_CCCC = 1200,
  NO_JET_MATCHING_SIG_COLLECTION = 1203,
  NO_JET_MATCHING_BKG_COLLECTION = 1204
};

class SaveFinalKinematics : public Processor, public TrueJet_Parser
{
	public:

		virtual Processor*  newProcessor()
		{
			return new SaveFinalKinematics;
		}
		SaveFinalKinematics();
		virtual ~SaveFinalKinematics() = default;
		SaveFinalKinematics(const SaveFinalKinematics&) = delete;
		SaveFinalKinematics& operator=(const SaveFinalKinematics&) = delete;
		virtual void init();
		virtual void Clear();
		virtual void processRunHeader( LCRunHeader*  /*run*/);
		virtual void processEvent( EVENT::LCEvent *pLCEvent );
		virtual void check();
		virtual void end();
		
 protected:
		void save_evt_with_error_code(int i_error_code);
		
		/**
		 * Add the expected output collections
		 */
		
		// Input collection names
		std::string m_inputMCTrueCollection{};
		std::string m_inputJetCollection{};
		std::string m_inputRecoParticleCollection{};
		std::string m_inputPreSelectionCollection{};

		// Inputs arguments
		std::string m_recoJetAlgoType{};
		std::string m_recoPfoAlgoType{};

		float m_minBLikeliness{};
		float m_minCLikeliness{};

		std::string m_outputFile{};
		std::string m_outputTree{};

		int require_presel_pass{};

		int debug_print{};

		// Runtime variables
		int n_run;
        int n_evt;
		int error_code{};
		int passed_preselection{};

		// List of IDs in MCParticlesSkimmed collection in following order: µ-, µ+, B1dec1, B1dec2, B2dec1, B2dec2 (where B is either Z or H) 
		const vector<int> mcp_ids_zhh = {8, 9, 12, 13, 14, 15};
		const vector<int> mcp_ids_zzh = {8, 9, 11, 10, 13, 14};

		vector<int> pdg_mcp{};
		vector<int> pdg_recopfo{};
		vector<int> pdg_recojet{};
		vector<int> pdg_truejet{};

		vector<vector<double>> fm_mcp{};
		vector<vector<double>> fm_recopfo{};
		vector<vector<double>> fm_recojet{};
		vector<vector<double>> fm_truejet{};

		vector<float> charge_mcp{};
		vector<float> charge_recopfo{};
		vector<float> charge_recojet{};
		vector<float> charge_truejet{};

		vector<int> truejet_jettypes{};

		TFile *m_pTFile{};
        TTree *m_pTTree{};
		
};

#endif
