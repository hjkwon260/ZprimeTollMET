#ifndef Trigger_H
#define Trigger_H

#include <vector>

namespace Ntuple
{
    class TriggerResult 
    {
      public:
        TriggerResult():
        isFired(false),
        Flag_goodVertices(false), Flag_globalSuperTightHalo2016Filter(false), Flag_HBHENoiseFilter(false), Flag_HBHENoiseIsoFilter(false), Flag_EcalDeadCellTriggerPrimitiveFilter(false),
        name("")
        {}
        ~TriggerResult(){}
        
        bool  isFired;
        bool  Flag_goodVertices, Flag_globalSuperTightHalo2016Filter, Flag_HBHENoiseFilter, Flag_HBHENoiseIsoFilter, Flag_EcalDeadCellTriggerPrimitiveFilter; //MET filter PAT
        std::string name; 
        
    };
    class TriggerObject 
    {
      public:
        TriggerObject():
        eta(0), phi(0),
        name("")
        {}
        ~TriggerObject(){}
        
        float eta, phi;
        std::string name;        
    };
    // class Event
    // {
    //     public:
    //     std::vector<Ntuple::GenParticle> genparticles_;    
    //     Event(){};
    //     ~Event(){};

    // };
    typedef std::vector<Ntuple::TriggerResult> TrigRes;
    typedef std::vector<Ntuple::TriggerObject> TrigObj;
}

// typedef std::vector<Ntuple::GenParticle> GenParticleCollection;

#endif
