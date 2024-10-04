#pragma once

#include "sndMuFilterBaseCut.h"

#include "TChain.h"

namespace snd{
  namespace analysis_cuts {
  
    class USPlanesHit : public MuFilterBaseCut {
    private :
      std::vector<int> planes_hit;
    public :
      USPlanesHit(std::vector<int> planes, TChain * tree);
      ~USPlanesHit(){;}
      bool passCut();
    };
    
  }
}
