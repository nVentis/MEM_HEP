#ifndef ReadAllKinematics_h
#define ReadAllKinematics_h 1

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

class ReadAllKinematics : public Processor, public TrueJet_Parser
{
	public:

		virtual Processor*  newProcessor()
		{
			return new ReadAllKinematics;
		}
		ReadAllKinematics();
		virtual ~ReadAllKinematics() = default;
		ReadAllKinematics(const ReadAllKinematics&) = delete;
		ReadAllKinematics& operator=(const ReadAllKinematics&) = delete;
		virtual void init();
		virtual void Clear();
		virtual void processRunHeader( LCRunHeader*  /*run*/);
		virtual void processEvent( EVENT::LCEvent *pLCEvent );
		virtual void check();
		virtual void end();
		
 protected:
		void save_evt_with_error_code(int error_code);
		
		/**
		 * Add the expected output collections
		 */
		
		/** Input collection name.
		 */
		std::string m_inputMCTrueCollection{};
		std::string m_inputLepPairCollection{};
		std::string m_inputJetCollection{};
		std::string m_inputPreSelectionCollection{};
		std::string m_inputHiggsPairCollection{};

		std::string m_inputJetMatchingSigCollection{};
		std::string m_inputJetMatchingBkgCollection{};

		std::string m_outputFile{};
		std::string m_outputTree{};

		LCCollection *inputJetCol{};

		int n_run;
        int n_evt;
		int error_code{};
		int passed_preselection{};
		int require_presel_pass{};

		// List of IDs in MCParticlesSkimmed collection in following order: µ-, µ+, B1dec1, B1dec2, B2dec1, B2dec2 (where B is either Z or H) 
		const vector<int> mcp_ids_zhh = {8, 9, 12, 13, 14, 15};
		const vector<int> mcp_ids_zzh = {8, 9, 11, 10, 13, 14};

		vector<vector<double>> fm_mcp{};
		vector<vector<double>> fm_recojet{};
		vector<vector<double>> fm_truejet{};

		vector<int> pdgs_mcp{};
		vector<int> pdgs_recojet{};
		vector<int> pdgs_truejet{};

		vector<int> recojet_pid_used{};
		vector<int> truejet_jettypes{};

		TFile *m_pTFile{};
        TTree *m_pTTree{};
		
};

#endif
