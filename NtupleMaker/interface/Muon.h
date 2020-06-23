#ifndef MUON_H
#define MUON_H

#include <vector>

namespace Ntuple
{
  class Muon 
  {
    public:
      Muon():
      isGlobalMuon(false), isStandAloneMuon(false), isTrackerMuon(false), isPFMuon(false),
      isTightMuon(false), isLooseMuon(false), isMediumMuon(false), isSoftMuon(false), isHighPtMuon(false), isNewHighPtMuon(false),
      pt(0), px(0), py(0), pz(0), eta(0), phi(0),
      charge(0), nChambers(0), stationMask(0), nMatchedStations(0),
      // normalizedChi2, nValidHits, nValidTrackerHits, nValidPixelHits, nTrackerLayers, qoverp, theta, 
      // staPt(0), staEta(0), staPhi(0),
      // pfPt(0), pfEta(0), pfPhi(0),
      trackIso(-1), ecalIso(-1), hcalIso(-1),
      pfchHadIso(-1), pfneuHadIso(-1), pfgammaIso(-1), pfpuIso(-1),
      // chHadIso03(-1), gammaIso03(-1), neuHadIso03(-1), puIso03(-1),
      // puppiChHadIso(-1), puppiGammaIso(-1), puppiNeuHadIso(-1), 
      // puppiChHadIsoNoLep(-1), puppiGammaIsoNoLep(-1), puppiNeuHadIsoNoLep(-1), 
      // d0(-999.), dz(-999.), sip3d(-999.),
      // tkNchi2(-999.), muNchi2(-999.),
      // trkKink(0), glbKink(0),
      // trkHitFrac(0), chi2LocPos(-999.), segComp(-999.), caloComp(-999.),
      // charge(0),
      // nValidHits(0),
      // typeBits(0), selectorBits(0), pogIDBits(0),
      // nTkHits(0), nPixHits(0),
      // nTkLayers(0), nPixLayers(0),
      // nMatchStn(0),
      // trkID(-1),
      // btt(-1),
      // hltMatchBits(0)
      //https://github.com/cms-sw/cmssw/blob/CMSSW_9_4_X/DataFormats/MuonReco/interface/Muon.h#L188-L212
      CutBasedIdLoose(false), CutBasedIdMedium(), CutBasedIdTight(), CutBasedIdGlobalHighPt(), CutBasedIdTrkHighPt(),
      PFIsoVeryLoose(), PFIsoLoose(), PFIsoMedium(), PFIsoTight(), PFIsoVeryTight(), TkIsoLoose(), TkIsoTight() 

      {}
      ~Muon(){}
      
      bool           isGlobalMuon, isStandAloneMuon, isTrackerMuon, isPFMuon;
      bool           isTightMuon, isLooseMuon, isMediumMuon, isSoftMuon, isHighPtMuon, isNewHighPtMuon;
      float          pt, px, py, pz, eta, phi;
      int            charge, nChambers, stationMask, nMatchedStations;       
      // float          staPt, staEta, staPhi;                 // STA track kinematics
      // float          pfPt, pfEta, pfPhi;                    // matched PFCandidate
      float          trackIso, ecalIso, hcalIso;              // detector isolation (R=0.3)
      float          pfchHadIso, pfneuHadIso, pfgammaIso, pfpuIso;  // PF isolation variables (R=0.4)
      // float          chHadIso03, gammaIso03, neuHadIso03, puIso03;  // PF isolation variables (R=0.3)
      // float          puppiChHadIso,      puppiGammaIso,      puppiNeuHadIso;  // Puppi Isolation R=0.4
      // float          puppiChHadIsoNoLep, puppiGammaIsoNoLep, puppiNeuHadIsoNoLep; // Puppi Isolation R=0.4 no lep
      // float          d0, dz, sip3d;                         // impact parameter
      // float          tkNchi2, muNchi2;                      // track fit normalized chi-square
      // float          trkKink, glbKink;                      // track kink
      // float          trkHitFrac;                            // fraction of valid tracker hits
      // float          chi2LocPos;                            // TRK-STA position match
      // float          segComp;                               // compatibility of tracker track with muon segment
      // float          caloComp;                              // muon hypothesis compatibility with calo energy
      // int            nValidHits;                            // number of valid muon hits in global fit
      // unsigned int   typeBits;                              // muon type bits
      // unsigned int   selectorBits;                          // MuonSelector bits
      // unsigned int   pogIDBits;                             // POG muon IDs from CMSSW
      // unsigned int   nTkHits, nPixHits;                     // number of hits in tracker
      // unsigned int   nTkLayers, nPixLayers;                 // number of hit layers in tracker
      // unsigned int   nMatchStn;                             // number of stations with muon segments
      // int            trkID;                                 // tracker track ID (unique per event)
      // int            btt;                                   // best track type
      // TriggerObjects hltMatchBits;                          // HLT matching
      bool CutBasedIdLoose, CutBasedIdMedium, CutBasedIdTight, CutBasedIdGlobalHighPt, CutBasedIdTrkHighPt;
      bool PFIsoVeryLoose, PFIsoLoose, PFIsoMedium, PFIsoTight, PFIsoVeryTight, TkIsoLoose, TkIsoTight;
          
  };

  typedef std::vector<Ntuple::Muon> MuonCollection;

  // enum EMuType
  // {
  //   // following convention in DataFormats/MuonReco/interface/Muon.h
  //   kGlobal     = 2,
  //   kTracker    = 4,
  //   kStandalone = 8,
  //   kCaloMuon   = 16,
  //   kPFMuon     = 32,
  //   kRPCMuon    = 64
  // };

  // enum EMuSelectorBit
  // { 
  //   // descriptions from DataFormats/MuonReco/interface/MuonSelectors.h
  //   kAll  			      = 0x0000001,  // dummy options - always true
  //   kAllGlobalMuons		      = 0x0000002,  // checks isGlobalMuon flag
  //   kAllStandAloneMuons		      = 0x0000004,  // checks isStandAloneMuon flag
  //   kAllTrackerMuons		      = 0x0000008,  // checks isTrackerMuon flag
  //   kTrackerMuonArbitrated	      = 0x0000010,  // resolve ambiguity of sharing segments
  //   kAllArbitrated		      = 0x0000020,  // all muons with the tracker muon arbitrated
  //   kGlobalMuonPromptTight	      = 0x0000040,  // global muons with tighter fit requirements
  //   kTMLastStationLoose		      = 0x0000080,  // penetration depth loose selector
  //   kTMLastStationTight		      = 0x0000100,  // penetration depth tight selector
  //   kTM2DCompatibilityLoose	      = 0x0000200,  // likelihood based loose selector
  //   kTM2DCompatibilityTight	      = 0x0000400,  // likelihood based tight selector
  //   kTMOneStationLoose		      = 0x0000800,  // require one well matched segment
  //   kTMOneStationTight		      = 0x0001000,  // require one well matched segment
  //   kTMLastStationOptimizedLowPtLoose = 0x0002000,  // combination of TMLastStation and TMOneStation
  //   kTMLastStationOptimizedLowPtTight = 0x0004000,  // combination of TMLastStation and TMOneStation
  //   kGMTkChiCompatibility 	      = 0x0008000,  // require tk stub have good chi2 relative to glb track
  //   kGMStaChiCompatibility	      = 0x0010000,  // require sta stub have good chi2 compatibility relative to glb track
  //   kGMTkKinkTight		      = 0x0020000,  // require a small kink value in the tracker stub
  //   kTMLastStationAngLoose	      = 0x0040000,  // TMLastStationLoose with additional angular cuts
  //   kTMLastStationAngTight	      = 0x0080000,  // TMLastStationTight with additional angular cuts
  //   kTMOneStationAngLoose 	      = 0x0100000,  // TMOneStationLoose with additional angular cuts
  //   kTMOneStationAngTight 	      = 0x0200000,  // TMOneStationTight with additional angular cuts
  //   //The two algorithms that follow are identical to what were known as
  //   //TMLastStationOptimizedLowPt* (sans the Barrel) as late as revision
  //   //1.7 of this file. The names were changed because indeed the low pt
  //   //optimization applies only to the barrel region, whereas the sel-
  //   //ectors above are more efficient at low pt in the endcaps, which is
  //   //what we feel is more suggestive of the algorithm name. This will be
  //   //less confusing for future generations of CMS members, I hope...
  //   kTMLastStationOptimizedBarrelLowPtLoose = 0x0400000,  // combination of TMLastStation and TMOneStation but with low pT optimization in barrel only
  //   kTMLastStationOptimizedBarrelLowPtTight = 0x0800000,  // combination of TMLastStation and TMOneStation but with low pT optimization in barrel only
  //   kRPCMuLoose                             = 0x1000000   // checks isRPCMuon flag (require two well matched hits in different RPC layers)
  // };

  // enum EPOGIDBit
  // {
    // Muon POG selection ID
  //   kPOGLooseMuon  =  1,
  //   kPOGMediumMuon =  2,
  //   kPOGTightMuon  =  4,
  //   kPOGSoftMuon   =  8,
  //   kPOGHighPtMuon = 16
  // };
}
#endif