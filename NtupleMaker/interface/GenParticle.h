#ifndef GENPARTICLE_H
#define GENPARTICLE_H

#include <vector>

namespace Ntuple
{
    class GenParticle 
    {
      public:
        GenParticle():
        mother(-1), pdgId(0),status(0), charge(0),
        pt(0), px(0), py(0), pz(0), eta(0), phi(0), energy(0), mass(0),
        isDirectHardProcessTauDecayProductFinalState(false), isDirectPromptTauDecayProductFinalState(false), isHardProcess(false), isLastCopy(false), isLastCopyBeforeFSR(false),
        fromHardProcessDecayed(false), fromHardProcessFinalState(false), fromHardProcessBeforeFSR(false), isPromptDecayed(false), isPromptFinalState(false)
        {}
        ~GenParticle(){}

        int   mother, pdgId, status, charge;
        float pt, px, py, pz, eta, phi, energy, mass;
        bool  isDirectHardProcessTauDecayProductFinalState, isDirectPromptTauDecayProductFinalState, isHardProcess, isLastCopy, isLastCopyBeforeFSR;
        bool  fromHardProcessDecayed, fromHardProcessFinalState, fromHardProcessBeforeFSR, isPromptDecayed, isPromptFinalState;
    };

    // class Event
    // {
    //     public:
    //     std::vector<Ntuple::GenParticle> genparticles_;    
    //     Event(){};
    //     ~Event(){};

    // };
    typedef std::vector<Ntuple::GenParticle> GenParticleCollection;
}

// typedef std::vector<Ntuple::GenParticle> GenParticleCollection;

#endif

