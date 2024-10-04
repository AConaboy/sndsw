#pragma once

#include "sndSciFiBaseCut.h"

#include "TChain.h"
#include "sndScifiHit.h"

#include "TVector3.h"
#include "TGraph.h"

namespace snd {
  namespace analysis_cuts {

    class scifiAngle : public sciFiBaseCut {
    private :
      double intercept_, slope_, max_chi2_;
      TVector3 a_, b_;
      TGraph * gv_;
      TGraph * gh_;
    public :
      scifiAngle(double intercept, double slope, double max_chi2, TChain * ch);
      ~scifiAngle(){;}

      bool passCut();

    };
  };
}
