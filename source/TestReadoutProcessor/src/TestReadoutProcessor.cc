#include "TestReadoutProcessor.h"
#include <iostream>

#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"

#ifdef MARLIN_USE_AIDA
#include <marlin/AIDAProcessor.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/ICloud1D.h>
//#include <AIDA/IHistogram1D.h>
#endif // MARLIN_USE_AIDA


using namespace lcio ;
using namespace marlin ;


TestReadoutProcessor aTestReadoutProcessor ;


TestReadoutProcessor::TestReadoutProcessor() : Processor("TestReadoutProcessor") {

    // modify processor description
    _description = "TestReadoutProcessor does whatever it does ..." ;


    // register steering parameters: name, description, class-variable, default value
    registerInputCollection( LCIO::MCPARTICLE,
            "CollectionName" , 
            "Name of the MCParticle collection"  ,
            _colName ,
            std::string("MCParticle")
    );
}



void TestReadoutProcessor::init() { 

    streamlog_out(DEBUG) << "   init called  " << std::endl ;

    // usually a good idea to
    printParameters() ;

    _nRun = 0 ;
    _nEvt = 0 ;

}


void TestReadoutProcessor::processRunHeader( LCRunHeader* run) { 

    _nRun++ ;
} 



void TestReadoutProcessor::processEvent( LCEvent * evt ) {
    
    //-- note: this will not be printed if compiled w/o MARLINDEBUG=1 !

    streamlog_out(DEBUG) << "   processing event: " << evt->getEventNumber() 
        << "   in run:  " << evt->getRunNumber() << std::endl ;



    _nEvt ++ ;
}



void TestReadoutProcessor::check( LCEvent * evt ) { 
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void TestReadoutProcessor::end(){ 

    //   std::cout << "TestReadoutProcessor::end()  " << name() 
    // 	    << " processed " << _nEvt << " events in " << _nRun << " runs "
    // 	    << std::endl ;

}