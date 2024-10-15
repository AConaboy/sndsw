#pragma once

#include "sndEventHeaderBaseCut.h"

#include "TChain.h"

namespace snd {
  namespace analysis_cuts {

    class ip1Cut : public snd::analysis_cuts::EventHeaderBaseCut {
    private:
    public :
      ip1Cut(TChain * ch);
      ~ip1Cut(){;}
      bool passCut();
    };

  }
}
