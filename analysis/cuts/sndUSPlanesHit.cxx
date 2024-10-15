#include "sndUSPlanesHit.h"

#include "TClonesArray.h"
#include "TChain.h"
#include "MuFilterHit.h"
#include "TString.h"

#include <vector>
#include <map>
#include <numeric>

namespace snd {
  namespace analysis_cuts {
  
    USPlanesHit::USPlanesHit(std::vector<int> planes, TChain * tree) : MuFilterBaseCut(tree) {
      planes_hit = planes;
      cutName = "Planes hit ";
      for (int plane : planes_hit) cutName += " "+std::to_string(plane);
      
      shortName = "USPlanesHit";
      nbins = std::vector<int>{1};
      range_start = std::vector<double>{0};
      range_end = std::vector<double>{1};
      plot_var = std::vector<double>{-1};
    }
    
    bool USPlanesHit::passCut(){
      MuFilterHit * hit;
      TIter hitIterator(muFilterDigiHitCollection);
      
      std::vector<bool> us = std::vector<bool>(planes_hit.size(), false);
      
      while ( (hit = (MuFilterHit*) hitIterator.Next()) ){
	if (! hit->isValid()) continue;

	if (hit->GetSystem() == 2) {
	  for (int i_plane = 0; i_plane < planes_hit.size(); i_plane++) {
	    if (hit->GetPlane() == planes_hit[i_plane]) us[i_plane] = true;
	  }
	}
      }
      
      for (bool this_plane_is_hit : us){
	if (not this_plane_is_hit){
	  return false;
	}
      }
      return true;
    }

  }
}
