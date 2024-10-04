#pragma once

#include <vector>

#include "sndBaseCut.h"

#include "TChain.h"
#include "TClonesArray.h"
#include "sndScifiHit.h"
#include "Scifi.h"
#include "SNDLHCEventHeader.h"
#include "TH1D.h"

namespace snd {
  namespace analysis_cuts {
  
    class sciFiBaseCut : public snd::analysis_cuts::baseCut {

    private : 
      static TChain * tree;
      static SNDLHCEventHeader * header;
      static unsigned long int read_entry;
      bool select_events_;
      double min_clock_cycle_;
      double max_clock_cycle_;
      double TDC2ns_;
      
    protected :
      static std::vector<TH1D *> hHitTime;

      static TClonesArray * scifiDigiHitCollection;
      static TClonesArray * scifiDigiHitCollection_raw;
      
      static std::vector<int> hits_per_plane_vertical;
      static std::vector<int> hits_per_plane_horizontal;

      static Scifi * scifiDet;
      
      void initializeEvent();

      sciFiBaseCut(TChain * ch, bool select_events = true,  double min_clock_cycle = -0.5, double max_clock_cycle = 1.2);
      ~sciFiBaseCut();
    };

  }
}
