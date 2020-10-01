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
        pt(0), px(0), py(0), phi(0), sumEt(0),
        puppi_pt(0), puppi_px(0), puppi_py(0), puppi_phi(0), puppi_sumEt(0),
        deep_pt(0), deep_px(0), deep_py(0), deep_phi(0), deep_sumEt(0)
        {}
        ~MET(){}
        
        float genpt, genpx, genpy, genphi, gensumEt;
        float pt, px, py, phi, sumEt;
        float puppi_pt, puppi_px, puppi_py, puppi_phi, puppi_sumEt;
        float deep_pt, deep_px, deep_py, deep_phi, deep_sumEt;
        
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

