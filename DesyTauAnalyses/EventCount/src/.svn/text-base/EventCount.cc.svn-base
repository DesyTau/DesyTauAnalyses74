// -*- C++ -*-
//
// Package:    EventCount
// Class:      EventCount
// 
/**\class EventCount EventCount.cc TrackAnalyzer/EventCount/src/EventCount.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Aliaksei Raspiareza,,,
//         Created:  Tue Feb  2 09:28:11 CET 2010
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <TH1.h>

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

//
// class declaration
//

class EventCount : public edm::EDAnalyzer {
   public:
      explicit EventCount(const edm::ParameterSet&);
      ~EventCount();


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      // ----------member data ---------------------------

  TH1F * _eventCount;
  TH1F * _eventCountJetBins;
  bool _doMC;
  edm::InputTag _genParticleSrc;


};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
EventCount::EventCount(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed

  _doMC            = iConfig.getParameter<bool>            (  "DoMC" );
  _genParticleSrc  = iConfig.getParameter<edm::InputTag>  (  "GenParticleSource" );


  edm::Service<TFileService> fs;

  _eventCount = fs->make<TH1F>("EventCount","number of events",1,-0.5,0.5);
  if (_doMC) 
    _eventCountJetBins = fs->make<TH1F>("EventCountJetBins","number of events in jet bins",11,-0.5,10.5);


}


EventCount::~EventCount()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
EventCount::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace reco;

   float dummy = 0.;
   _eventCount->Fill(dummy);

   if (_doMC) {
     Handle<GenParticleCollection> genParticles;
     iEvent.getByLabel(_genParticleSrc,genParticles);
     int nParticles = genParticles->size();
     // Counting jets --->
     int Npartons = 0;
     bool count_jets = false;
     for (int iP=0;iP<nParticles;++iP) {
       GenParticleRef part(genParticles,iP);
       if (part->status() != 3) continue;
       int pdg = abs(part->pdgId());
       if (count_jets) {
	 if (pdg == 1 || pdg == 2 || pdg == 3 || pdg == 4 || pdg == 5 || pdg == 6 || pdg == 21) Npartons++;
       }
       if (pdg == 23) count_jets = true;  // Start counting partons after we find the Z

     }
     // end of counting jets
     _eventCountJetBins->Fill(float(Npartons));
   }


#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void 
EventCount::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
EventCount::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(EventCount);
