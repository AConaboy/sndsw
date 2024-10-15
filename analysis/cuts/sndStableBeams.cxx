#include "sndStableBeams.h"

#include "SNDLHCEventHeaderConst.h"

#include <iostream>

namespace snd::analysis_cuts {

  stableBeamsCut::stableBeamsCut(TChain * ch) : EventHeaderBaseCut(ch) {
    cutName = "Stable beams";

    shortName = "StableBeams";
    nbins = std::vector<int>{31};
    range_start = std::vector<double>{0};
    range_end = std::vector<double>{31};
    plot_var = std::vector<double>{-1};

  }

  bool stableBeamsCut::passCut(){

    plot_var[0] = header->GetBeamMode();
    
    if (header->GetBeamMode() == (int) LhcBeamMode::StableBeams) return true;
    else return false;
  }
}
