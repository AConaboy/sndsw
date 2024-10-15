#include "sndSciFiContinuity.h"

#include "sndSciFiTools.h"

#include "TChain.h"

namespace snd {
  namespace analysis_cuts {

    sciFiContinuity::sciFiContinuity(TChain * ch) : sciFiBaseCut(ch){
      cutName = "SciFi continuity";

      shortName = "SciFiContinuity";
      nbins = std::vector<int>{1};
      range_start = std::vector<double>{0};
      range_end = std::vector<double>{1};
      plot_var = std::vector<double>{-1};

    }

    bool sciFiContinuity::passCut(){
      initializeEvent();
    
      bool gotFirstPlane = false;
      // For a plane to count, need both planes to have hits
      for (int i_plane = 0; i_plane < hits_per_plane_horizontal.size(); i_plane++){
	if ((hits_per_plane_horizontal[i_plane]*hits_per_plane_vertical[i_plane] == 0) and gotFirstPlane){
	  return false;
	} else if (hits_per_plane_horizontal[i_plane]*hits_per_plane_vertical[i_plane] > 0){
	  gotFirstPlane = true;
	}
      }
      return true;
    }
  }	     
}
