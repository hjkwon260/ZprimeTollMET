// -*- C++ -*-
//
// Package:    ZprimeTollMET/NtupleMaker
// Class:      NtupleMaker
// 
/**\class NtupleMaker NtupleMaker.cc ZprimeTollMET/NtupleMaker/plugins/NtupleMaker.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Hyejin Kwon
//         Created:  Sat, 09 Nov 2019 08:43:08 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
//////////////////////////
// -- Track & Vertex -- //
//////////////////////////
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "PhysicsTools/RecoUtils/interface/CandCommonVertexFitter.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
////////////////////
// -- Triggers -- //
////////////////////
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
////////////////////////////
// -- For GenParticles -- //
////////////////////////////
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
/////////////////////
// -- For Muons -- //
/////////////////////
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"

/////////////////////////
// -- For Electrons -- //
/////////////////////////
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "RecoEgamma/EgammaTools/interface/EffectiveAreas.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
////////////////////
// -- For Jets -- //
////////////////////
#include "DataFormats/PatCandidates/interface/Jet.h"
///////////////////
// -- For MET -- //
///////////////////
#include "DataFormats/PatCandidates/interface/MET.h"

#include "ZprimeTollMET/NtupleMaker/interface/Trigger.h"
#include "ZprimeTollMET/NtupleMaker/interface/GenParticle.h"
#include "ZprimeTollMET/NtupleMaker/interface/Muon.h"
#include "ZprimeTollMET/NtupleMaker/interface/Electron.h"
#include "ZprimeTollMET/NtupleMaker/interface/Jet.h"
#include "ZprimeTollMET/NtupleMaker/interface/MET.h"
#include "ZprimeTollMET/NtupleMaker/interface/Event.h"

#include "TTree.h"
#include <vector>

using namespace std;
using namespace reco;
using namespace edm;

//
// class declaration
//


class NtupleMaker : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit NtupleMaker(const edm::ParameterSet&);
      ~NtupleMaker();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      virtual void fillEvent(const edm::Event& iEvent, bool isMC);
      virtual void fillTriggers(const edm::Event& iEvent);
      virtual void fillGenParticles(const edm::Event& iEvent);
      virtual void fillMuons(const edm::Event& iEvent);
      virtual void fillElectrons(const edm::Event& iEvent);
      virtual void fillMET(const edm::Event& iEvent);
      virtual void fillJet(const edm::Event& iEvent);
      bool isNewHighPtMuon(const reco::Muon& muon, const reco::Vertex& vtx);

      // ----------member data ---------------------------

      bool isMC;

      edm::EDGetTokenT<std::vector<PileupSummaryInfo>> t_pu;
      edm::EDGetTokenT<edm::TriggerResults> t_trigresult;
      edm::EDGetTokenT<edm::TriggerResults> t_trigresultRECO;
      edm::EDGetTokenT<edm::TriggerResults> t_trigresultPAT;
      // edm::EDGetTokenT<std::vector<pat::TriggerObjectStandAlone>> t_trigobject; 
      edm::EDGetTokenT<double> t_prefweight;
      edm::EDGetTokenT<double> t_prefweightup;
      edm::EDGetTokenT<double> t_prefweightdown;     
      edm::EDGetTokenT<std::vector<reco::Vertex>> t_pv;
      edm::EDGetTokenT<reco::GenParticleCollection> t_genparticle;
      edm::EDGetTokenT<GenEventInfoProduct> t_genevtinfo;
      edm::EDGetTokenT<std::vector<pat::Muon>> t_muon;
      edm::EDGetTokenT<std::vector<pat::Electron>> t_electron;
      edm::EDGetTokenT<std::vector<pat::Jet>> t_jet;
      // edm::EDGetTokenT<std::vector<pat::Jet>> t_DeepFlavour;
      edm::EDGetTokenT<std::vector<pat::MET>> t_MET;
      edm::EDGetTokenT<std::vector<pat::MET>> t_puppiMET;
      edm::EDGetTokenT<std::vector<pat::MET>> t_deepMET;
      
      // edm::Handle<std::vector<reco::Vertex>> h_pv;

      TTree *ntuple_;
      // reco::Vertex pv;
      Ntuple::EventInfo eventinfo_;
      Ntuple::TrigRes trigresults_;
      // Ntuple::TrigObj trigobjects_;
      Ntuple::GenParticleCollection genparticles_;
      Ntuple::MuonCollection muons_;
      Ntuple::ElectronCollection electrons_;
      Ntuple::JetCollection jets_;
      Ntuple::METevt metevt_;

      // double prefiringweight_ = 0;
      // double prefiringweightup_ = 0;
      // double prefiringweightdown_ = 0;
      int ngenparticles_ = 0;
      float genweight_ = 1;
      int nmuons_ = 0;
      int nelectrons_ = 0;
      int njets_ = 0;
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
NtupleMaker::NtupleMaker(const edm::ParameterSet& iConfig):
isMC(iConfig.getParameter<bool>("isMC")),
t_pu                (consumes<std::vector< PileupSummaryInfo >>          (iConfig.getUntrackedParameter<edm::InputTag>("pileupsummary"))),
t_trigresult        (consumes<edm::TriggerResults>                       (iConfig.getUntrackedParameter<edm::InputTag>("triggerresults"))),
t_trigresultRECO     (consumes<edm::TriggerResults>                      (iConfig.getUntrackedParameter<edm::InputTag>("triggerresultsRECO"))),
t_trigresultPAT     (consumes<edm::TriggerResults>                       (iConfig.getUntrackedParameter<edm::InputTag>("triggerresultsPAT"))),
// t_trigobject        (consumes<std::vector<pat::TriggerObjectStandAlone>> (iConfig.getUntrackedParameter<edm::InputTag>("triggerobjects"))),
t_prefweight        (consumes<double>                                    (iConfig.getUntrackedParameter<edm::InputTag>("prefiringweight"))),
t_prefweightup      (consumes<double>                                    (iConfig.getUntrackedParameter<edm::InputTag>("prefiringweightUp"))),
t_prefweightdown    (consumes<double>                                    (iConfig.getUntrackedParameter<edm::InputTag>("prefiringweightDown"))),
t_pv                (consumes<std::vector<reco::Vertex>>                 (iConfig.getUntrackedParameter<edm::InputTag>("pvs"))),
t_genparticle       (consumes<reco::GenParticleCollection>               (iConfig.getUntrackedParameter<edm::InputTag>("genparticles"))),
t_genevtinfo        (consumes<GenEventInfoProduct>                       (iConfig.getUntrackedParameter<edm::InputTag>("genevtinfo"))),
t_muon              (consumes<std::vector<pat::Muon>>                    (iConfig.getUntrackedParameter<edm::InputTag>("muons"))),
t_electron          (consumes<std::vector<pat::Electron>>                (iConfig.getUntrackedParameter<edm::InputTag>("electrons"))),
t_jet               (consumes<std::vector<pat::Jet>>                     (iConfig.getUntrackedParameter<edm::InputTag>("jets"))),
t_MET               (consumes<std::vector<pat::MET>>                     (iConfig.getUntrackedParameter<edm::InputTag>("MET"))),
t_puppiMET          (consumes<std::vector<pat::MET>>                     (iConfig.getUntrackedParameter<edm::InputTag>("puppiMET"))),
t_deepMET           (consumes<std::vector<pat::MET>>                     (iConfig.getUntrackedParameter<edm::InputTag>("deepMET")))

{}


NtupleMaker::~NtupleMaker()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void NtupleMaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  // iEvent.getByToken(t_pv, h_pv);
  // pv = h_pv->front();

  eventinfo_.clear();
  trigresults_.clear();
  // trigobjects_.clear();
  genparticles_.clear();
  muons_.clear();
  electrons_.clear();
  jets_.clear();
  metevt_.clear();
  // event_.genparticles_.clear();

  fillEvent(iEvent, isMC);
  fillTriggers(iEvent);
  if(isMC) fillGenParticles(iEvent);
  fillMuons(iEvent);
  fillElectrons(iEvent);
  // fillJet(iEvent);
  // fillMET(iEvent);
  ntuple_->Fill();
 
  // ntuple_->Branch("genparticles" ,&event_.genparticles_);
}


// ------------ method called once each job just before starting event loop  ------------
void NtupleMaker::beginJob()
{
  edm::Service<TFileService> fs;
  ntuple_ = fs->make<TTree>("ntuple","ntuple");  

  ntuple_->Branch("eventinfo" ,&eventinfo_);
  ntuple_->Branch("trigresults" ,&trigresults_);
  // ntuple_->Branch("trigobjects" ,&trigobjects_);
  // ntuple_->Branch("prefiringweight" ,&prefiringweight_, "prefiringweight/D");
  // ntuple_->Branch("prefiringweightup" ,&prefiringweightup_, "prefiringweightup/D");
  // ntuple_->Branch("prefiringweightdown" ,&prefiringweightdown_, "prefiringweightdown/D");
  ntuple_->Branch("ngenparticles" ,&ngenparticles_, "ngenparticles/I");
  ntuple_->Branch("genweight" ,&genweight_, "genweight/F");
  ntuple_->Branch("genparticles" ,&genparticles_);
  ntuple_->Branch("nmuons" ,&nmuons_, "nmuons/I");
  ntuple_->Branch("muons" ,&muons_);  
  ntuple_->Branch("nelectrons" ,&nelectrons_, "nelectrons/I");
  ntuple_->Branch("electrons" ,&electrons_); 
  ntuple_->Branch("njets" ,&njets_, "njets/I");  
  ntuple_->Branch("jets" ,&jets_);  
  ntuple_->Branch("met" ,&metevt_); 
}

// ------------ method called once each job just after ending the event loop  ------------
void NtupleMaker::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
NtupleMaker::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

void NtupleMaker::fillEvent(const edm::Event &iEvent, bool isMC) {

  Ntuple::Event evt_;

  evt_.runNum       = iEvent.id().run();
  evt_.lumi = iEvent.id().luminosityBlock();
  evt_.eventNum     = iEvent.id().event();

  edm::Handle<std::vector<reco::Vertex>> h_pv;
  iEvent.getByToken(t_pv, h_pv);

  evt_.nvertices = h_pv->size();

  if(isMC){
    edm::Handle<std::vector<PileupSummaryInfo>>  h_pu;
    iEvent.getByToken(t_pu, h_pu);

    for(std::vector<PileupSummaryInfo>::const_iterator pu = h_pu->begin(); pu != h_pu->end(); ++pu)
    {
      int BX = pu->getBunchCrossing();

      if( BX == 0 )
      {
        // npvin = PVI->getPU_NumInteractions(); // in time only
        evt_.npileup = pu->getTrueNumInteractions(); // in and out of time
        continue;
      }
    }
  }  

  edm::Handle<double> h_prefweight;
  iEvent.getByToken(t_prefweight, h_prefweight);
  evt_.prefiringweight = (*h_prefweight);

  edm::Handle<double> h_prefweightup;
  iEvent.getByToken(t_prefweightup, h_prefweightup);
  evt_.prefiringweightup = (*h_prefweightup);

  edm::Handle<double> h_prefweightdown;
  iEvent.getByToken(t_prefweightdown, h_prefweightdown);
  evt_.prefiringweightdown = (*h_prefweightdown);

  eventinfo_.push_back(evt_);
  
}

void NtupleMaker::fillTriggers(const edm::Event &iEvent) {

    Ntuple::TriggerResult trigres_;
    // Ntuple::TriggerObject trigobj_;

    edm::Handle< edm::TriggerResults > h_trigresults;
    edm::Handle< edm::TriggerResults > h_trigresultsPAT;
    iEvent.getByToken(t_trigresult,h_trigresults);
    // if(isMC) iEvent.getByToken(t_trigresultPAT,h_trigresultsPAT);
    iEvent.getByToken(t_trigresultPAT,h_trigresultsPAT);
    if(!(h_trigresultsPAT.isValid())) iEvent.getByToken(t_trigresultRECO,h_trigresultsPAT);

    // edm::Handle<std::vector<pat::TriggerObjectStandAlone>> h_trigobjects;
    // iEvent.getByToken(t_trigobject, h_trigobjects);

    string trigs[] = {
        // "HLT_IsoMu17_eta2p1_v*",
        // "HLT_IsoMu20_v*",
        // "HLT_IsoTkMu20_v*", 
        "HLT_IsoMu24_v*",
        "HLT_IsoTkMu24_v*",
    //     "HLT_IsoTkMu20_eta2p1_v*", 
    // "HLT_DoubleEle33_CaloIdL_MW_v*",
    // "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v*",
    // "HLT_DoubleEle37_Ele27_CaloIdL_GsfTrkIdVL_v*",
        "HLT_Mu50_v*",
    "HLT_TkMu50_v*",
    "HLT_Ele27_WPTight_Gsf_v*",
    "HLT_Ele32_WPTight_Gsf_v*",
    "HLT_Ele32_WPMedium_Gsf_v*",
    //     "HLT_Mu45_eta2p1_v*",       
    //     "HLT_Mu17_TrkIsoVVL_v*",
    //     "HLT_DoubleIsoMu17_eta2p1_v*",                
    //     "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*",
    //     "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*",
    //     "HLT_Mu27_TkMu8_v*",
    // "HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_v*",
    // "HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*",
    // "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v*",
    // "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*",
    //     "HLT_Ele22_eta2p1_WP75_Gsf_v*",
    //     "HLT_Ele22_eta2p1_WPLoose_Gsf_v*",
    //     "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v*",
    // "Mu"
    };

    const unsigned ntrigs = sizeof(trigs)/sizeof(*trigs);
    const unsigned ntrigres = h_trigresults.product()->size();

    edm::TriggerNames triglist = iEvent.triggerNames(*h_trigresults);

    // cout << h_trigresultsPAT.product()->size( << endl;

    //test filter
    // if(isMC){
    edm::TriggerNames triglistPAT = iEvent.triggerNames(*h_trigresultsPAT);
    trigres_.Flag_goodVertices = h_trigresultsPAT.product()->accept(triglistPAT.triggerIndex("Flag_goodVertices"));    
    trigres_.Flag_globalSuperTightHalo2016Filter = h_trigresultsPAT.product()->accept(triglistPAT.triggerIndex("Flag_globalSuperTightHalo2016Filter"));    
    trigres_.Flag_HBHENoiseFilter = h_trigresultsPAT.product()->accept(triglistPAT.triggerIndex("Flag_HBHENoiseFilter"));    
    trigres_.Flag_HBHENoiseIsoFilter = h_trigresultsPAT.product()->accept(triglistPAT.triggerIndex("Flag_HBHENoiseIsoFilter"));    
    trigres_.Flag_EcalDeadCellTriggerPrimitiveFilter = h_trigresultsPAT.product()->accept(triglistPAT.triggerIndex("Flag_EcalDeadCellTriggerPrimitiveFilter"));  
    // }  

    for( unsigned i = 0; i != ntrigres; i++ ) {
        std::string _trigname = triglist.triggerName(i);

        // cout << _trigname << endl;

        for( unsigned j = 0; j != ntrigs; j++ ) {
            if( _trigname.find(trigs[j].substr(0,trigs[j].find("*"))) != std::string::npos ) {
                trigres_.name = _trigname;

                trigres_.isFired = (h_trigresults.product()->accept(i) ? true : false);
                trigresults_.push_back(trigres_);
            }
        }
    }

    // for(unsigned i=0; i<h_trigobjects->size(); i++) {

    //     pat::TriggerObjectStandAlone trigobj = h_trigobjects->at(i);
    //     trigobj.unpackPathNames(triglist);
    //     const std::vector< std::string > pathnames = trigobj.pathNames();

    //     for(std::vector<std::string>::const_iterator name = pathnames.begin(); name != pathnames.end(); ++name){
    //         for(unsigned int j=0; j<ntrigs; j++) {
    //             if( name->find(trigs[j].substr(0,trigs[j].find("*"))) != std::string::npos && trigobj.hasPathName(*name,true,true) ) {
    //                 trigobj_.eta = trigobj.eta();
    //                 trigobj_.phi = trigobj.phi();
    //                 trigobj_.name = *name;
    //                 trigobjects_.push_back(trigobj_);
    //             }
    //         }
    //     }
    // }

}

void NtupleMaker::fillGenParticles(const edm::Event& iEvent) {

    Ntuple::GenParticle genpar_;

    edm::Handle<reco::GenParticleCollection> h_genparticles;
    iEvent.getByToken(t_genparticle, h_genparticles);

    // edm::Handle<GenEventInfoProduct> genInfo;
    // iEvent.getByToken(generatorToken, genInfo);

    // if( genInfo->weight() > 0 ) event_.weight = 1.0;
    // else event_.weight = -1.0;

    ngenparticles_ = h_genparticles->size();
    if( h_genparticles->size() == 0 ) return;

    int _ngenparticle = 0;
    for( reco::GenParticleCollection::const_iterator genpar = h_genparticles->begin(); genpar != h_genparticles->end(); ++genpar ) {

        if( abs(genpar->pdgId()) == 5 || abs(genpar->pdgId()) == 6 || abs(genpar->pdgId()) == 11 || abs(genpar->pdgId()) == 12 || abs(genpar->pdgId()) == 13 || abs(genpar->pdgId()) == 14 || abs(genpar->pdgId()) == 15 || abs(genpar->pdgId()) == 16 || abs(genpar->pdgId()) == 23 || abs(genpar->pdgId()) == 24 || abs(genpar->pdgId()) == 32 || abs(genpar->pdgId()) == 1000022 || abs(genpar->pdgId()) == 1000024 ) {

            genpar_.pdgId  = genpar->pdgId(); 
            genpar_.pt     = genpar->pt(); 
            genpar_.px     = genpar->px();
            genpar_.py     = genpar->py();
            genpar_.pz     = genpar->pz();
            genpar_.eta    = genpar->eta();
            genpar_.phi    = genpar->phi();
            genpar_.energy = genpar->energy();
            genpar_.mass   = genpar->mass();
            genpar_.charge = genpar->charge();
            genpar_.status = genpar->status();
            genpar_.mother = genpar->mother(0)->pdgId();

            genpar_.isDirectHardProcessTauDecayProductFinalState = genpar->isDirectHardProcessTauDecayProductFinalState();
            genpar_.isDirectPromptTauDecayProductFinalState      = genpar->isDirectPromptTauDecayProductFinalState();
            genpar_.isHardProcess                                = genpar->isHardProcess();
            genpar_.isLastCopy                                   = genpar->isLastCopy();
            genpar_.isLastCopyBeforeFSR                          = genpar->isLastCopyBeforeFSR();
            genpar_.isPromptDecayed                              = genpar->isPromptDecayed();
            genpar_.isPromptFinalState                           = genpar->isPromptFinalState();
            genpar_.fromHardProcessFinalState                    = genpar->fromHardProcessFinalState();
            genpar_.fromHardProcessDecayed                       = genpar->fromHardProcessDecayed();
            genpar_.fromHardProcessBeforeFSR                     = genpar->fromHardProcessBeforeFSR();

            _ngenparticle++;
            genparticles_.push_back(genpar_);
            // event_.genparticles_.push_back(genpar_);
        }
    }
    ngenparticles_ = _ngenparticle;

    edm::Handle<GenEventInfoProduct> h_genevtinfo;
    iEvent.getByToken(t_genevtinfo, h_genevtinfo);
    genweight_ = h_genevtinfo->weight();
    // if(h_genevtinfo->weight()<0) genweight_=-1;
    // cout << genweight_ << endl;
}

void NtupleMaker::fillMuons(const edm::Event& iEvent) {

    Ntuple::Muon mu_;

    edm::Handle<std::vector<reco::Vertex>> h_pv;
    iEvent.getByToken(t_pv, h_pv);

    reco::Vertex pv = h_pv->front();

    edm::Handle<std::vector<pat::Muon>> h_muons;
    iEvent.getByToken(t_muon,h_muons);

    // edm::Handle<pat::PackedCandidateCollection> pfcands;
    // iEvent.getByToken(pfCandidatesToken, pfcands);

    // edm::Handle< double > miniRhoH;
    // iEvent.getByToken(miniRhoToken,miniRhoH);

    nmuons_ = h_muons->size();
    if( h_muons->size() == 0 ) return;

    // int _nmuons = 0;
    for( std::vector<pat::Muon>::const_iterator mu = h_muons->begin(); mu != h_muons->end(); ++mu ) {    

        // const pat::Muon &lep = *mu;

        mu_.isStandAloneMuon = mu->isStandAloneMuon();
        mu_.isGlobalMuon     = mu->isGlobalMuon();
        mu_.isTrackerMuon    = mu->isTrackerMuon();
        mu_.isPFMuon         = mu->isPFMuon();
        mu_.isTightMuon      = mu->isTightMuon(pv);
        mu_.isLooseMuon      = mu->isLooseMuon();
        mu_.isMediumMuon     = mu->isMediumMuon(); 
        mu_.isSoftMuon       = mu->isSoftMuon(pv);
        mu_.isHighPtMuon     = mu->isHighPtMuon(pv); 
        mu_.px               = mu->px();
        mu_.py               = mu->py();
        mu_.pz               = mu->pz();
        mu_.pt               = mu->pt();
        mu_.eta              = mu->eta();
        mu_.phi              = mu->phi();
        mu_.charge           = mu->charge();
        mu_.nChambers        = mu->numberOfChambers();
        mu_.stationMask      = mu->stationMask();
        mu_.nMatchedStations = mu->numberOfMatchedStations();
        if( mu->innerTrack().isNonnull() ) mu_.nTkLayers  = mu->innerTrack()->hitPattern().trackerLayersWithMeasurement();
        mu_.isNewHighPtMuon = isNewHighPtMuon( (*mu), pv );

        // if( mu->isGlobalMuon() ) { // Global Muon

        //     reco::TrackRef glbTrack  =  mu->globalTrack();

        //     if( glbTrack.isNonnull() ) {

        //         const reco::HitPattern & glbhit = glbTrack->hitPattern();

        //         mu_.normalizedChi2  =  glbTrack->normalizedChi2();
        //         mu_.nValidHits      =  glbTrack->numberOfValidHits();         
        //         mu_.nValidMuonHits  = glbhit.numberOfValidMuonHits();
        //         mu_.qoverp          = glbTrack->qoverp();
        //         mu_.theta           = glbTrack->theta();
        //         mu_.lambda          = glbTrack->lambda();
        //         mu_.dxy             = glbTrack->dxy();
        //         mu_.d0              = glbTrack->d0();
        //         mu_.dsz             = glbTrack->dsz();
        //         mu_.dz              = glbTrack->dz();
        //         mu_.dxyBS           = glbTrack->dxy(beamSpot.position());
        //         mu_.dszBS           = glbTrack->dsz(beamSpot.position());
        //         mu_.dzBS            = glbTrack->dz(beamSpot.position());
        //         mu_.vx              = glbTrack->vx();
        //         mu_.vy              = glbTrack->vy();
        //         mu_.vz              = glbTrack->vz();

        //     }

        //     reco::TrackRef trackerTrack = mu->innerTrack();

        //     if( trackerTrack.isNonnull() ) {

        //         const reco::HitPattern & inhit = trackerTrack->hitPattern();

        //         mu_.nValidTrackerHits = inhit.numberOfValidTrackerHits();
        //         mu_.nValidPixelHits   = inhit.numberOfValidPixelHits();
        //         mu_.nTrackerLayers    = inhit.trackerLayersWithMeasurement();

        //     }

        // } // Global Muon 
        // else if( mu->isStandAloneMuon() ) { // STA Muon

        //     reco::TrackRef muonTrack  =  mu->outerTrack();

        //     if( muonTrack.isNonnull() ) {

        //         const reco::HitPattern & muonhit  =  muonTrack->hitPattern();

        //         mu_.nValidMuonHits  =  muonhit.numberOfValidMuonHits();

        //     }

        // } // STA Muon
        // else if(mu->isTrackerMuon() ) { // Tracker Muon

        //     reco::TrackRef trackerTrack = mu->innerTrack();

        //     if( trackerTrack.isNonnull() ) {

        //         const reco::HitPattern & inhit = trackerTrack->hitPattern();

        //         mu_.normalizedChi2    = trackerTrack->normalizedChi2();
        //         mu_.nValidHits        = trackerTrack->numberOfValidHits();                    
        //         mu_.nValidTrackerHits = inhit.numberOfValidTrackerHits();
        //         mu_.nValidPixelHits   = inhit.numberOfValidPixelHits();
        //         mu_.nTrackerLayers    = inhit.trackerLayersWithMeasurement();
        //         mu_.qoverp            = trackerTrack->qoverp();
        //         mu_.theta             = trackerTrack->theta();
        //         mu_.lambda            = trackerTrack->lambda();
        //         mu_.dxy               = trackerTrack->dxy();
        //         mu_.d0                = trackerTrack->d0();
        //         mu_.dsz               = trackerTrack->dsz();
        //         mu_.dz                = trackerTrack->dz();
        //         mu_.dxyBS             = trackerTrack->dxy(beamSpot.position());
        //         mu_.dszBS             = trackerTrack->dsz(beamSpot.position());
        //         mu_.dzBS              = trackerTrack->dz(beamSpot.position());
        //         mu_.vx                = trackerTrack->vx();
        //         mu_.vy                = trackerTrack->vy();
        //         mu_.vz                = trackerTrack->vz();

        //     }
        // } // Tracker Muon

        // // MuonBestTrack
        // if( mu->muonBestTrack().isNonnull() ) {

        //     mu_.muonBestTrack_nTrackerLayers = mu->muonBestTrack()->hitPattern().trackerLayersWithMeasurement();
        //     mu_.muonBestTrack_px             = mu->muonBestTrack()->px();
        //     mu_.muonBestTrack_py             = mu->muonBestTrack()->py();
        //     mu_.muonBestTrack_pz             = mu->muonBestTrack()->pz();
        //     mu_.muonBestTrack_pt             = mu->muonBestTrack()->pt();
        //     mu_.muonBestTrack_ptError        = mu->muonBestTrack()->ptError();

        //     if( !pvHandle->empty() && !pvHandle->front().isFake() ) {

        //         mu_.dxyVTX = mu->muonBestTrack()->dxy(vtx.position());
        //         mu_.dszVTX = mu->muonBestTrack()->dsz(vtx.position());
        //         mu_.dzVTX  = mu->muonBestTrack()->dz(vtx.position());
        //     }

        // }

        // // InnerTrack
        // if( mu->innerTrack().isNonnull() ) {

        //     mu_.innerTrack_nTrackerLayers = mu->innerTrack()->hitPattern().trackerLayersWithMeasurement();
        //     mu_.innerTrack_px             = mu->innerTrack()->px();
        //     mu_.innerTrack_py             = mu->innerTrack()->py();
        //     mu_.innerTrack_pz             = mu->innerTrack()->pz();
        //     mu_.innerTrack_pt             = mu->innerTrack()->pt();
        //     mu_.innerTrack_ptError        = mu->innerTrack()->ptError();

        //     if( !pvHandle->empty() && !pvHandle->front().isFake() ) {

        //         mu_.innerTrack_dxyVTX = mu->innerTrack()->dxy(vtx.position());
        //         mu_.innerTrack_dszVTX = mu->innerTrack()->dsz(vtx.position());
        //         mu_.innerTrack_dzVTX  = mu->innerTrack()->dz(vtx.position());

        //     }

        // }

        // // TunePMuonBestTrack
        // if( mu->tunePMuonBestTrack().isNonnull() ) {

        //     mu_.tunePMuonBestTrack_nTrackerLayers = mu->tunePMuonBestTrack()->hitPattern().trackerLayersWithMeasurement();
        //     mu_.tunePMuonBestTrack_px             = mu->tunePMuonBestTrack()->px();
        //     mu_.tunePMuonBestTrack_py             = mu->tunePMuonBestTrack()->py();
        //     mu_.tunePMuonBestTrack_pz             = mu->tunePMuonBestTrack()->pz();
        //     mu_.tunePMuonBestTrack_pt             = mu->tunePMuonBestTrack()->pt();
        //     mu_.tunePMuonBestTrack_ptError        = mu->tunePMuonBestTrack()->ptError();

        //     if( !pvHandle->empty() && !pvHandle->front().isFake() ) {

        //         mu_.tunePMuonBestTrack_dxyVTX = mu->tunePMuonBestTrack()->dxy(vtx.position());
        //         mu_.tunePMuonBestTrack_dszVTX = mu->tunePMuonBestTrack()->dsz(vtx.position());
        //         mu_.tunePMuonBestTrack_dzVTX  = mu->tunePMuonBestTrack()->dz(vtx.position());

        //     }

        // }   

        // Isolation
        mu_.trackIso    = mu->trackIso();
        mu_.ecalIso     = mu->ecalIso();
        mu_.hcalIso     = mu->hcalIso();
        mu_.pfchHadIso  = mu->pfIsolationR04().sumChargedHadronPt;
        mu_.pfneuHadIso = mu->pfIsolationR04().sumNeutralHadronEt;
        mu_.pfgammaIso  = mu->pfIsolationR04().sumPhotonEt;
        mu_.pfpuIso     = mu->pfIsolationR04().sumPUPt;

        //from selector
        mu_.CutBasedIdLoose        = mu->passed(pat::Muon::CutBasedIdLoose);
        mu_.CutBasedIdMedium       = mu->passed(pat::Muon::CutBasedIdMedium); 
        mu_.CutBasedIdTight        = mu->passed(pat::Muon::CutBasedIdTight);
        mu_.CutBasedIdGlobalHighPt = mu->passed(pat::Muon::CutBasedIdGlobalHighPt);
        mu_.CutBasedIdTrkHighPt    = mu->passed(pat::Muon::CutBasedIdTrkHighPt);
        mu_.PFIsoVeryLoose         = mu->passed(pat::Muon::PFIsoVeryLoose);
        mu_.PFIsoLoose             = mu->passed(pat::Muon::PFIsoLoose);
        mu_.PFIsoMedium            = mu->passed(pat::Muon::PFIsoMedium);
        mu_.PFIsoTight             = mu->passed(pat::Muon::PFIsoTight);
        mu_.PFIsoVeryTight         = mu->passed(pat::Muon::PFIsoVeryTight);
        mu_.TkIsoLoose             = mu->passed(pat::Muon::TkIsoLoose);
        mu_.TkIsoTight             = mu->passed(pat::Muon::TkIsoTight);

        // _nmuons++;
        muons_.push_back(mu_);
    }
    // nmuons_ = _nmuons;
}

void NtupleMaker::fillElectrons(const edm::Event& iEvent) {

    Ntuple::Electron el_;

    edm::Handle<std::vector<pat::Electron>> h_electrons;
    iEvent.getByToken(t_electron,h_electrons);

    nelectrons_ = h_electrons->size();
    if( h_electrons->size() == 0 ) return;

    for( std::vector<pat::Electron>::const_iterator el = h_electrons->begin(); el != h_electrons->end(); ++el ) {    

        el_.CutBasedIdVeto   = el->electronID("cutBasedElectronID-Fall17-94X-V2-veto");
        el_.CutBasedIdLoose  = el->electronID("cutBasedElectronID-Fall17-94X-V2-loose");
        el_.CutBasedIdMedium = el->electronID("cutBasedElectronID-Fall17-94X-V2-medium");
        el_.CutBasedIdTight  = el->electronID("cutBasedElectronID-Fall17-94X-V2-tight");
        el_.MVAwp80iso       = el->electronID("mvaEleID-Fall17-iso-V2-wp80");
        el_.MVAwp80noiso     = el->electronID("mvaEleID-Fall17-noIso-V2-wp80");
        el_.MVAwp90iso       = el->electronID("mvaEleID-Fall17-iso-V2-wp90");
        el_.MVAwp90noiso     = el->electronID("mvaEleID-Fall17-noIso-V2-wp90");
        el_.uncorrpx               = el->px();
        el_.uncorrpy               = el->py();
        el_.uncorrpz               = el->pz();
        el_.uncorrpt               = el->pt();
        el_.uncorreta              = el->eta();
        el_.uncorrphi              = el->phi();
        el_.uncorrenergy           = el->energy();
        el_.charge           = el->charge();
        // auto corrP4  = pat::Electron::p4() * pat::Electron::userFloat("ecalTrkEnergyPostCorr") / pat::Electron::energy();
        auto corrp4 = el->p4()*el->userFloat("ecalTrkEnergyPostCorr")/el->energy();
        el_.px           = corrp4.px();
        el_.py           = corrp4.py();
        el_.pz           = corrp4.pz();
        el_.pt           = corrp4.pt();
        el_.eta          = corrp4.eta();
        el_.phi          = corrp4.phi();
        el_.energy       = corrp4.energy();       

        // Isolation
        el_.pfchHadIso  = el->pfIsolationVariables().sumChargedHadronPt;
        el_.pfneuHadIso = el->pfIsolationVariables().sumNeutralHadronEt;
        el_.pfgammaIso  = el->pfIsolationVariables().sumPhotonEt;
        el_.pfpuIso     = el->pfIsolationVariables().sumPUPt;
        el_.pfrelIso    = (el_.pfchHadIso + max<float>( 0.0, el_.pfneuHadIso + el_.pfgammaIso - 0.5 * el_.pfpuIso))/(el->pt());

        electrons_.push_back(el_);
    }
}

void NtupleMaker::fillJet(const edm::Event& iEvent) {

    Ntuple::Jet jet_;

    edm::Handle<std::vector<pat::Jet>> h_jets;
    iEvent.getByToken(t_jet,h_jets);

    njets_ = h_jets->size();

    if( h_jets->size() == 0 ) return;

    for( vector<pat::Jet>::const_iterator jet = h_jets->begin(); jet != h_jets->end(); ++jet ) {

        jet_.pt       = jet->pt();
        jet_.eta      = jet->eta();
        jet_.phi      = jet->phi();
        jet_.charge   = jet->jetCharge();
        jet_.flavor   = jet->partonFlavour();
        jet_.btag     = jet->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
        jet_.pfDeepCSV_probb     = jet->bDiscriminator("pfDeepCSVJetTags:probb");
        jet_.pfDeepCSV_probbb     = jet->bDiscriminator("pfDeepCSVJetTags:probbb");
        jet_.pfDeepCSV     = jet->bDiscriminator("pfDeepCSVJetTags:probbb")+jet->bDiscriminator("pfDeepCSVJetTags:probb");
        jet_.pfDeepFlv     = jet->bDiscriminator("pfDeepFlavourJetTags:probb")+jet->bDiscriminator("pfDeepFlavourJetTags:probbb")+jet->bDiscriminator("pfDeepFlavourJetTags:problepb");
        jet_.chfrac   = jet->chargedHadronEnergyFraction();
        jet_.nhfrac   = jet->neutralHadronEnergyFraction();
        jet_.nhemfrac = jet->neutralEmEnergyFraction();
        jet_.chemfrac = jet->chargedEmEnergyFraction();
        jet_.chmulti  = jet->chargedMultiplicity();
        jet_.nhmulti  = jet->neutralMultiplicity();

        jets_.push_back(jet_);
    }

}

void NtupleMaker::fillMET(const edm::Event& iEvent) {

    Ntuple::MET met_;

    edm::Handle<std::vector<pat::MET>> met;
    iEvent.getByToken(t_MET,met);

    if(isMC)
    {
      met_.genpt          = met->front().genMET()->pt();
      met_.genphi         = met->front().genMET()->phi();
      met_.genpx          = met->front().genMET()->px();
      met_.genpy          = met->front().genMET()->py();
      met_.gensumEt       = met->front().genMET()->sumEt();
    }
    met_.pt    = met->front().pt();
    met_.phi   = met->front().phi();
    met_.px    = met->front().px();
    met_.py    = met->front().py();
    met_.sumEt = met->front().sumEt();

    edm::Handle<std::vector<pat::MET>> puppimet;
    iEvent.getByToken(t_puppiMET,puppimet);

    met_.puppi_pt    = puppimet->front().pt();
    met_.puppi_phi   = puppimet->front().phi();
    met_.puppi_px    = puppimet->front().px();
    met_.puppi_py    = puppimet->front().py();
    met_.puppi_sumEt = puppimet->front().sumEt();

    edm::Handle<std::vector<pat::MET>> deepmet;
    iEvent.getByToken(t_deepMET,deepmet);

    met_.deep_pt    = deepmet->front().pt();
    met_.deep_phi   = deepmet->front().phi();
    met_.deep_px    = deepmet->front().px();
    met_.deep_py    = deepmet->front().py();
    met_.deep_sumEt = deepmet->front().sumEt();

    metevt_.push_back(met_);

}

//until CMSSW_10_4_X
bool NtupleMaker::isNewHighPtMuon(const reco::Muon& muon, const reco::Vertex& vtx){
  if(!muon.isGlobalMuon()) return false;

  bool muValHits = ( muon.globalTrack()->hitPattern().numberOfValidMuonHits()>0 ||
                     muon.tunePMuonBestTrack()->hitPattern().numberOfValidMuonHits()>0 );

  bool muMatchedSt = muon.numberOfMatchedStations()>1;
  if(!muMatchedSt) {
    if( muon.isTrackerMuon() && muon.numberOfMatchedStations()==1 ) {
      if( muon.expectedNnumberOfMatchedStations()<2 ||
          !(muon.stationMask()==1 || muon.stationMask()==16) ||
          muon.numberOfMatchedRPCLayers()>2
        )
        muMatchedSt = true;
    }
  }

  bool muID = muValHits && muMatchedSt;

  bool hits = muon.innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5 &&
    muon.innerTrack()->hitPattern().numberOfValidPixelHits() > 0;

  bool momQuality = muon.tunePMuonBestTrack()->ptError()/muon.tunePMuonBestTrack()->pt() < 0.3;

  bool ip = fabs(muon.innerTrack()->dxy(vtx.position())) < 0.2 && fabs(muon.innerTrack()->dz(vtx.position())) < 0.5;

  return muID && hits && momQuality && ip;

}
//define this as a plug-in
DEFINE_FWK_MODULE(NtupleMaker);
