/**
 * @author Bryan Bliewert (TUM/DESY) <bryan.bliewert@nventis.eu>
 * @date 13.12.2023
 *
 */

#ifndef JetConvProcessor_h
#define JetConvProcessor_h 1

#include "torch/torch.h"
#include "torch/script.h"

#include <pybind11/embed.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <memory>
#include <iostream>

#include "marlin/Processor.h"
#include <EVENT/LCCollection.h>
#include "lcio.h"
#include <string>
#include <vector>
#include "TLorentzVector.h"
#include "TrueJet_Parser.h"

using namespace lcio;
using namespace marlin;
using namespace pybind11::literals;
namespace py = pybind11;

class JetConvProcessor : public Processor, public TrueJet_Parser
{
public:
	virtual Processor *newProcessor()
	{
		return new JetConvProcessor;
	}
	JetConvProcessor();
	virtual ~JetConvProcessor() = default;
	JetConvProcessor(const JetConvProcessor &) = delete;
	JetConvProcessor &operator=(const JetConvProcessor &) = delete;
	virtual void init();
	virtual void clear();
	virtual void processRunHeader(LCRunHeader * /*run*/);
	virtual void processEvent(EVENT::LCEvent *pLCEvent);
	virtual void check();
	virtual void end();

protected:
	/**
	 * Input arguments
	 * */
	std::string m_inputMCTrueCollection{};
	std::string m_inputPfoCollection{};
	int m_nJets{};
	std::string m_torchScriptPath{};
	std::string m_outputJetCollection{};

	/**
	 * 
	 * 
	*/
	std::string m_inputLepPairCollection{};
	std::string m_inputJetCollection{};
	std::string m_inputHdecayModeCollection{};
	std::string m_inputHiggsPairCollection{};

	std::string m_inputJetMatchingSigCollection{};
	std::string m_inputJetMatchingBkgCollection{};

	int m_nRun;
	int m_nEvt;
	int m_error_code;

	// Torch
	torch::jit::script::Module m_model{};
	torch::TensorOptions m_options{};

	// pybind
	py::object m_spectral_clustering{};
};

#endif
