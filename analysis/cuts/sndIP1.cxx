#include "sndIP1.h"

#include "SNDLHCEventHeaderConst.h"

#include <iostream>

namespace snd::analysis_cuts {

  ip1Cut::ip1Cut(TChain * ch) : EventHeaderBaseCut(ch) {
    cutName = "IP1 bunch crossing";

    shortName = "IP1";

    /* Plot variable encoded like this:
       0 isIP1();
       1 isB1Only();
       2 isB2noB1();
    */
    nbins = std::vector<int>{3};
    range_start = std::vector<double>{0};
    range_end = std::vector<double>{3};
    plot_var = std::vector<double>{-1};

  }

  bool ip1Cut::passCut(){

    bool ret = false;
    if (header->isIP1()){
      ret = true;
      plot_var[0] = 0;
    } else if (header->isB1Only()){
      plot_var[0] = 1;
    } else if (header->isB2noB1()){
      plot_var[0] = 2;
    }
    return ret;
  }
}
