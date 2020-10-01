#ifndef Event_H
#define Event_H

#include <vector>

namespace Ntuple
{
    class Event 
    {
      public:
        Event():
        nvertices(-1), npileup(-1),
        prefiringweight(0), prefiringweightup(0), prefiringweightdown(0)

        {}
        ~Event(){}
        
        int nvertices, npileup;
        double prefiringweight, prefiringweightup, prefiringweightdown;
        
    };

    typedef std::vector<Ntuple::Event> EventInfo;
}

#endif

