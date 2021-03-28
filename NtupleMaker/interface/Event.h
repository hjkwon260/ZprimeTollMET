#ifndef Event_H
#define Event_H

#include <vector>

namespace Ntuple
{
    class Event 
    {
      public:
        Event():
        runNum(-1), lumi(-1), eventNum(0),
        nvertices(-1), npileup(-1),
        prefiringweight(0), prefiringweightup(0), prefiringweightdown(0)

        {}
        ~Event(){}
        
        int runNum, lumi;
        unsigned long long eventNum;
        int nvertices, npileup;
        double prefiringweight, prefiringweightup, prefiringweightdown;
        
    };

    typedef std::vector<Ntuple::Event> EventInfo;
}

#endif

