#ifndef MET_H
#define MET_H

#include <vector>

namespace Ntuple
{
    class MET 
    {
      public:
        MET():
        genpt(0), genpx(0), genpy(0), genphi(0), gensumEt(0),
        pt(0), px(0), py(0), phi(0), sumEt(0)
        {}
        ~MET(){}
        
        float genpt, genpx, genpy, genphi, gensumEt;
        float pt, px, py, phi, sumEt;
        
    };

    // class Event
    // {
    //     public:
    //     std::vector<Ntuple::GenParticle> genparticles_;    
    //     Event(){};
    //     ~Event(){};

    // };
    typedef std::vector<Ntuple::MET> METevt;
}

// typedef std::vector<Ntuple::GenParticle> GenParticleCollection;

#endif

