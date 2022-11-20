#ifndef Jet_H
#define Jet_H

#include <vector>

namespace Ntuple
{
    class Jet 
    {
      public:
        Jet():
        pt(0), eta(0), phi(0),
        charge(0), flavor(0),
        btag(0), pfDeepCSV_probb(0), pfDeepCSV_probbb(0), pfDeepCSV(0), pfDeepFlv(0),
        chfrac(0), nhfrac(0), nhemfrac(0), chemfrac(0), 
        chmulti(0), nhmulti(0)
        {}
        ~Jet(){}
        
        float pt, eta, phi;
        int charge, flavor;
        float btag, pfDeepCSV_probb, pfDeepCSV_probbb, pfDeepCSV, pfDeepFlv;
        float chfrac, nhfrac,nhemfrac, chemfrac;
        int chmulti, nhmulti;
        
    };

    // class Event
    // {
    //     public:
    //     std::vector<Ntuple::GenParticle> genparticles_;    
    //     Event(){};
    //     ~Event(){};

    // };
    typedef std::vector<Ntuple::Jet> JetCollection;
}

// typedef std::vector<Ntuple::GenParticle> GenParticleCollection;

#endif

