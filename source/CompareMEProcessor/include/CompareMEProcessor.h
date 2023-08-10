/**
 * Used to calculate MEs for ZHH and ZZH processes based on pure MC truth data, without TrueJet-jet clustering.
 * TODO: Include flag for using TrueJet output
 *
 * @author Bryan Bliewert (TUM/DESY) <bryan.bliewert@tum.de>
 * @date 20.05.2023
 *
 */

#ifndef CompareMEProcessor_h
#define CompareMEProcessor_h 1

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
#include "TrueJet_Parser.h"
class TFile;
class TH1F;
class TH1I;
class TH2I;
class TTree;

using namespace lcio;
using namespace marlin;
using namespace lcme;

enum ERRORS : unsigned int
{
	PRESELECTION_FAILED_BUT_REQUIRED = 1,
	UNHANDLED_PROCESS = 5,
	TRUEJET_NULL = 7,
	TRUEJET_NO_JETS = 8,

	INCOMPLETE_RECO_LEPTON_PAIR = 101,
	INCOMPLETE_RECOJET_COLLECTION = 102,

	INCOMPLETE_TRUEJET_LEPTON_PAIR = 201,
	INCOMPLETE_TRUEJET_COLLECTION = 202,

	SOME_TRUEJET_NOT_FOUND = 306,



	NEITHER_BBBB_NOR_CCCC = 1200,
	NO_JET_MATCHING_SIG_COLLECTION = 1203,
	NO_JET_MATCHING_BKG_COLLECTION = 1204
};

class CompareMEProcessor : public Processor, public TrueJet_Parser
{
public:
	virtual Processor *newProcessor()
	{
		return new CompareMEProcessor;
	}
	CompareMEProcessor();
	virtual ~CompareMEProcessor() = default;
	CompareMEProcessor(const CompareMEProcessor &) = delete;
	CompareMEProcessor &operator=(const CompareMEProcessor &) = delete;
	virtual void init();
	virtual void Clear();
	virtual void processRunHeader(LCRunHeader * /*run*/);
	virtual void processEvent(EVENT::LCEvent *pLCEvent);
	virtual void check();
	virtual void end();
	int getZDecayModeFromPDG(int pdg);

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
	std::string m_inputHdecayModeCollection{};
	std::string m_inputHiggsPairCollection{};

	std::string m_inputJetMatchingSigCollection{};
	std::string m_inputJetMatchingBkgCollection{};

	std::string m_outputFile{};
	std::string m_outputTree{};

	int m_nRun;
	int m_nEvt;
	int m_mode_me;
	int m_error_code{};

	int m_mode{}; // 0 => use MCTruth data; 1 => use reconstructed data (HiggsPair, LeptonPair, HdecayMode, and some jet e.g. RefinedJets)
	int m_lepton_mode{};
	int m_truejet_mode{}; // 0 => True; 1 => Seen; (2 for True-of-Seen probably not required)
	int m_saveTransferEnergies{};
	int m_saveTransferKinematics{};
	int m_require_presel_pass{};
	float m_Hmass{};

	TFile *m_pTFile{};
	TTree *m_pTTree{};
	lcme::LCMEZHH *_zhh; // ZHH MEM calculator instance
	lcme::LCMEZZH *_zzh; // ZZH MEM calculator instance

	// Control parameters
	int m_zhh_is_set{}; // used internally
	int m_zzh_is_set{}; // used internally

	// True
	int m_true_h1_decay_pdg{}; // abs(PDG) of particle H1 decayed to
	int m_true_h2_decay_pdg{}; // abs(PDG) of particle H2 decayed to (!= 0 only for true ZHH events)
	int m_true_z2_decay_pdg{}; // abs(PDG) of particle Z2 decayed to (!= 0 only for true ZZH events)

	// partons_{1,2,3,4} only != 0 if there indeed were quarks
	int m_parton1_pdg{};
	int m_parton2_pdg{};
	int m_parton3_pdg{};
	int m_parton4_pdg{};
	
	float m_parton1_e{};
	float m_parton2_e{};
	float m_parton3_e{};
	float m_parton4_e{};

	std::vector<float> m_parton1_p{};

	float m_true_lep1_e{};
	float m_true_lep2_e{};

	// Event data
	int m_z1_decay_pdg{};  // as input from parameter; defaults to 5 (mu+mu-)
	int m_z1_decay_mode{}; // m_z1_decay_pdg converted to corresponding LCME value
	int m_z2_decay_mode{}; // when ZZH is assumed; from m_z2_decay1_pdg
	int m_is_zhh{};		   // true label
	int m_is_zzh{};		   // true label
	int m_passed_preselection{};

	int m_h1_decay_pdg{}; // abs(PDG) of particle H1 decayed to
	int m_h2_decay_pdg{}; // abs(PDG) of particle H2 decayed to (>0 only for true ZHH events)
	int m_z2_decay_pdg{}; // abs(PDG) of particle Z2 decayed to (>0 only for true ZZH events)

	// Only for reco and TrueJet data; thus only filled for m_mode=1
	float m_jet1_e{};
	float m_jet2_e{};
	float m_jet3_e{};
	float m_jet4_e{};

	float m_lep1_e{};
	float m_lep2_e{};

	// Di-jet properties
	int m_misclustering_region{}; //e.g. AA, BB; see Misclustering code
	int m_misclustering_region_icns{}; // same as for m_misclustering_region

	float m_efrac1_reco{};
	float m_efrac2_reco{};
	float m_efrac1_true{};
	float m_efrac2_true{};

	float m_efrac1_icn_reco{};
	float m_efrac2_icn_reco{};
	float m_efrac1_icn_true{};
	float m_efrac2_icn_true{};

	LCCollection *inputJetCol{};
	LCCollection *inputHiggsPair{};

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
	float m_zhh_q2_z{};
	float m_zhh_q2_h1{};
	float m_zhh_q2_h2{};

	// 1. Assuming ZZH
	// 1.a ZZH output
	float m_zzh_sigma{};
	float m_zzh_sigmalll{};
	float m_zzh_sigmallz{};
	float m_zzh_sigmallr{};

	float m_zzh_sigmalrl{};
	float m_zzh_sigmalrz{};
	float m_zzh_sigmalrr{};

	float m_zzh_sigmarrr{};
	float m_zzh_sigmarrz{};
	float m_zzh_sigmarrl{};

	float m_zzh_sigmarlr{};
	float m_zzh_sigmarlz{};
	float m_zzh_sigmarll{};

	float m_zzh_mz1{};
	float m_zzh_mz2{};
	float m_zzh_mzz{};
	float m_zzh_mzzh{};
	float m_zzh_mh{};

	float m_zzh_phi{};
	float m_zzh_phiz{};
	float m_zzh_phiz1f{};
	float m_zzh_phiz2f{};

	float m_zzh_costheta{};
	float m_zzh_costhetaz{};
	float m_zzh_costhetaz1f{};
	float m_zzh_costhetaz2f{};

	// 1.b ZZH input
	float m_zzh_q2_z1{};
	float m_zzh_q2_z2{};
	float m_zzh_q2_h{};
};

#endif
