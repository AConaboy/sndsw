#pragma once

#include "sndEventHeaderBaseCut.h"

#include "TChain.h"

namespace snd {
  namespace analysis_cuts {

    class stableBeamsCut : public snd::analysis_cuts::EventHeaderBaseCut {
    private:
    public :
      stableBeamsCut(TChain * ch);
      ~stableBeamsCut(){;}
      bool passCut();
    };

  }
}
