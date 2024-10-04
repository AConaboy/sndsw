#pragma once

#include "sndSciFiBaseCut.h"

#include "TChain.h"
#include "sndScifiHit.h"

namespace snd{
  namespace analysis_cuts{
    
    class sciFiContinuity : public sciFiBaseCut {
    public :
      sciFiContinuity(TChain * ch);
      ~sciFiContinuity(){;}

      bool passCut();

    };
  };
}
