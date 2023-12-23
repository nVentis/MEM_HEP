/**
 * @author Bryan Bliewert (TUM/DESY) <bryan.bliewert@nventis.eu>
 * @date 13.12.2023
*/

#include "JetConvProcessor.h"
#include <vector>
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>
#include <IMPL/ReconstructedParticleImpl.h>
#include <UTIL/PIDHandler.h>

using namespace lcio ;
using namespace marlin ;
namespace py = pybind11 ;

JetConvProcessor aJetConvProcessor ;

JetConvProcessor::JetConvProcessor() :

  Processor("JetConvProcessor"),
  m_nJets(0),
  m_nRun(0),
  m_nEvt(0),
  m_error_code(0)
{

	_description = "JetConvProcessor writes relevant observables to root-file " ;

  registerInputCollection(LCIO::MCPARTICLE,
				 "InputMCTrueCollection",
				 "all MCParticles collection",
				 m_inputMCTrueCollection,
				 std::string("MCParticlesSkimmed")
				 );

  registerInputCollection(LCIO::RECONSTRUCTEDPARTICLE,
				 "inputPfoCollection",
				 "all PandoraPFOs",
				 m_inputPfoCollection,
				 std::string("PandoraPFOs")
				 );

  registerProcessorParameter("nJets",
        "in how many jets to cluster",
        m_nJets,
        int(4)
        );

  registerProcessorParameter("torchScriptPath",
        "name of output root file",
        m_torchScriptPath,
        std::string("/nfs/dust/ilc/user/bliewert/jet_training/model_3m/best.pt")
        );

  registerOutputCollection(LCIO::RECONSTRUCTEDPARTICLE,
				 "outputJetCollection",
				 "output jet collection",
				 m_outputJetCollection,
				 std::string("GraphJets")
				 );
}

void JetConvProcessor::init()
{  
  streamlog_out(DEBUG) << "   init called  " << std::endl;this->clear();

  // Load TorchScript model
  streamlog_out(MESSAGE) << "Loading torch script model from \'" << m_torchScriptPath << "\'" << std::endl;

  try {   
    m_model = torch::jit::load(m_torchScriptPath);
    m_model.to(torch::kCPU);
    m_model.eval();
    // streamlog_out(MESSAGE) << m_model.dump_to_str(true, true, true) << '\n'; // print model parameters

    m_options = torch::TensorOptions().dtype(torch::kFloat32).device(torch::kCPU);

    streamlog_out(DEBUG) << "Successfully imported PyTorch model" << std::endl;
  } catch (const c10::Error &e) {
    std::stringstream message;
    message << "Error loading the PyTorch model from \'" << m_torchScriptPath << "\': " << e.what();
    streamlog_out(ERROR) << message.str() << std::endl;
    throw EVENT::Exception(message.str());
  }

  streamlog_out(MESSAGE) << "Loading python and sklearn" << std::endl;

  try {
    py::initialize_interpreter();
    py::module_ sklearn_cluster = py::module_::import("sklearn.cluster");
    m_spectral_clustering = sklearn_cluster.attr("SpectralClustering");

    streamlog_out(DEBUG) << "Successfully imported sklearn" << std::endl;
  } catch (int e) {
    streamlog_out(DEBUG) << "Error while trying to import sklearn" << std::endl;      
  }

  streamlog_out(DEBUG) << "init finished  " << std::endl;
}

void JetConvProcessor::clear() 
{
  streamlog_out(DEBUG) << "clear called  " << std::endl;
}

void JetConvProcessor::processRunHeader( LCRunHeader*  /*run*/) { 
  m_nRun++ ;
}

void JetConvProcessor::processEvent( EVENT::LCEvent *pLCEvent )
{
  this->clear();
  
  m_nRun = pLCEvent->getRunNumber();
  m_nEvt = pLCEvent->getEventNumber();
  streamlog_out(DEBUG) << "processing event: " << pLCEvent->getEventNumber() << " in run: " << pLCEvent->getRunNumber() << std::endl;

  try {
    // Fetching collections
    LCCollection *inputMCTrueCollection;
    LCCollection *inputPFOCollection;
    LCCollectionVec *outputJetCollection = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);

    streamlog_out(DEBUG) << "  getting true MC collection: " << m_inputMCTrueCollection << std::endl ;
    inputMCTrueCollection = pLCEvent->getCollection( m_inputMCTrueCollection );

    streamlog_out(DEBUG) << "  getting preselection_passed collection: " << m_inputPfoCollection << std::endl ;
    inputPFOCollection = pLCEvent->getCollection( m_inputPfoCollection );

    // Disable auto-grad
    torch::NoGradGuard no_grad;
    
    int n_pfos = inputPFOCollection->getNumberOfElements();
    auto input_pfos = torch::zeros({n_pfos, 4},
      torch::TensorOptions()
      .dtype(torch::kFloat32)
      .layout(torch::kStrided)
      .device(torch::kCPU));
      
    for (int i = 0; i < n_pfos; ++i) {
      const EVENT::ReconstructedParticle *pfo = dynamic_cast<EVENT::ReconstructedParticle*>(inputPFOCollection->getElementAt(i));
      auto momentum = pfo->getMomentum();
      input_pfos[i][0] = pfo->getEnergy();
      input_pfos[i][1] = momentum[0];
      input_pfos[i][2] = momentum[1];
      input_pfos[i][3] = momentum[2];
    }

    // Convert resuling tensor to C array
    constexpr size_t elsize = sizeof(float);
    torch::Tensor similarity_tensor = m_model.run_method("get_similarity", input_pfos).toTensor();
    //float* similarity_arr = similarity_tensor.data_ptr<float>();
    py::array_t<float> similarity_python(
      {n_pfos, n_pfos}, // shape
      {n_pfos*elsize, elsize}, // strides
      similarity_tensor.data_ptr<float>());

    streamlog_out(DEBUG) << "Calculated similarity matrix" << std::endl;

    // Debug print similarity
  #ifdef DEBUG
    std::cout << "Event " << m_nEvt << " | Found " << i << " PFOs" << std::endl;
    std::cout << "Shape " << input_pfos.size(0) << ":" << input_pfos.size(1) << std::endl;
    for (int j = 0; j < n_pfos; j++) {
      std::cout << similarity_tensor[j][0].item<float>() << "|";
      //std::cout << (*similarity_arr++) << "|";
    }
    std::cout << std::endl;
  #endif

    // Construct clusters
    py::dict kwargs = py::dict(
      "n_clusters"_a = m_nJets,
      "affinity"_a = "precomputed"
    );
    py::object SpectralClustering = m_spectral_clustering(**kwargs);

    streamlog_out(DEBUG) << "Initialized clustering" << std::endl;

    py::list cluster_id_vec_python = SpectralClustering.attr("fit_predict")(similarity_python);
    std::vector<int> cluster_id_vec = cluster_id_vec_python.cast<std::vector<int>>();

    streamlog_out(DEBUG) << "Calculated clusters" << std::endl;

  #ifdef DEBUG
    int k = 0;
    for (py::handle item: cluster_id_vec_python){
      std::cout << k++ << " : " << item.attr("__str__")().cast<std::string>() << std::endl;
    }
    std::cout << std::endl;
  #endif

    // Save to LCIO
    for (int j = 0; j < m_nJets; j++) {
      ReconstructedParticleImpl* graphjet = new ReconstructedParticleImpl;

      for (int i = 0; i < n_pfos; i++) {
        if (cluster_id_vec[i] == j) {
          EVENT::ReconstructedParticle* pfo = dynamic_cast<EVENT::ReconstructedParticle*>(inputPFOCollection->getElementAt(i));
          graphjet->addParticle(pfo);

          const double* jet_mom = graphjet->getMomentum();
          const double* pfo_mom = pfo->getMomentum();
          float new_mom[3];
          new_mom[0] = jet_mom[0] + pfo_mom[0];
          new_mom[1] = jet_mom[1] + pfo_mom[1];
          new_mom[2] = jet_mom[2] + pfo_mom[2];

          graphjet->setCharge(graphjet->getCharge() + pfo->getCharge());
          graphjet->setMomentum(new_mom);
          graphjet->setEnergy(graphjet->getEnergy() + pfo->getEnergy());
        }        
      }

      graphjet->setMass(std::sqrt(
        std::pow(graphjet->getEnergy(), 2) - (
          std::pow(graphjet->getMomentum()[0], 2) +
          std::pow(graphjet->getMomentum()[1], 2) +
          std::pow(graphjet->getMomentum()[2], 2)
        )
      ));
      graphjet->setType(4); // 0: unknown 1: single 2:v0 3: compound 4:jet
      //graphjet->setReferencePoint(vpos);

      outputJetCollection->addElement(graphjet);

      streamlog_out(DEBUG) << "Jet " << j << " -> " << (graphjet->getParticles()).size() << " PFOs" << std::endl;
    }

    pLCEvent->addCollection(outputJetCollection, m_outputJetCollection.c_str());

  } catch(DataNotAvailableException &e) {
    streamlog_out(MESSAGE) << "processEvent : Input collections not found in event " << m_nEvt << std::endl;
  }
}

void JetConvProcessor::check()
{
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}

void JetConvProcessor::end()
{
  py::finalize_interpreter();
  //this->delall();
}