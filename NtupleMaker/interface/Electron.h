#ifndef Electron_H
#define Electron_H

#include <vector>

namespace Ntuple
{
  class Electron 
  {
    public:
      Electron():
      CutBasedIdVeto(false), CutBasedIdLoose(false), CutBasedIdMedium(false), CutBasedIdTight(false),
      MVAwp80iso(false), MVAwp80noiso(false), MVAwp90iso(false), MVAwp90noiso(false),
      pt(0), px(0), py(0), pz(0), eta(0), phi(0), energy(0), 
      charge(0),
      pfchHadIso(-1), pfneuHadIso(-1), pfgammaIso(-1), pfpuIso(-1), pfrelIso(-1)

      {}
      ~Electron(){}
      
      bool           CutBasedIdVeto, CutBasedIdLoose, CutBasedIdMedium, CutBasedIdTight;
      bool           MVAwp80iso, MVAwp80noiso, MVAwp90iso, MVAwp90noiso;
      float          pt, px, py, pz, eta, phi, energy;
      int            charge;       
      float          pfchHadIso, pfneuHadIso, pfgammaIso, pfpuIso, pfrelIso;    // detector isolation (R=0.3)
          
  };

  typedef std::vector<Ntuple::Electron> ElectronCollection;


}
#endif