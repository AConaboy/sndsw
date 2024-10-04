#include "sndUSBarsVeto.h"

#include "TClonesArray.h"
#include "TChain.h"
#include "MuFilterHit.h"
#include "TString.h"

#include <vector>
#include <map>
#include <numeric>

namespace snd{
  namespace analysis_cuts {

    USBarsVeto::USBarsVeto(std::vector<std::pair<int, double> > avg_per_plane, bool bottom, TChain * tree) : MuFilterBaseCut(tree) {
      bars = avg_per_plane;
      bottom_ = bottom;
      cutName = "Bars to exclude ";
      for (std::pair<int, double> bar : bars) cutName += " ("+std::to_string(std::get<0>(bar))+" "+std::to_string(std::get<1>(bar))+")";

      shortName = "USBarsVeto";
      for (std::pair<int, double> bar : bars) shortName += "_"+std::to_string(std::get<0>(bar))+"_"+std::to_string(std::get<1>(bar));
      
      nbins = std::vector<int>(bars.size(), 20);
      range_start = std::vector<double>(bars.size(), 0);
      range_end = std::vector<double>(bars.size(), 10);
      plot_var = std::vector<double>(bars.size(), 0.);
    }

    bool USBarsVeto::passCut(){
      MuFilterHit * hit;
      TIter hitIterator(muFilterDigiHitCollection);

      std::fill(plot_var.begin(), plot_var.end(), 0.);
      std::vector<int> bars_hit(bars.size(), 0);
      
      while ( (hit = (MuFilterHit*) hitIterator.Next()) ){
	if (! hit->isValid()) continue;
      
	if (hit->GetSystem() == 2) {
	  int this_plane = hit->GetPlane();
	  int this_bar = int(hit->GetDetectorID()%100);

	  for (int i_plane = 0; i_plane < bars.size(); i_plane++){
	    if (this_plane == std::get<0>(bars.at(i_plane))) {
	      bars_hit.at(i_plane)++;
	      plot_var.at(i_plane) += this_bar;
	    }
	  }
	}
      }

      bool ret = true;
      for (int i_plane = 0; i_plane < bars.size(); i_plane++){
	plot_var.at(i_plane) /= bars_hit.at(i_plane);
	if (bottom_){	
	  if (plot_var.at(i_plane) < std::get<1>(bars.at(i_plane))) {
	    ret = false;
	  }
	} else {
	  if (plot_var.at(i_plane) >= std::get<1>(bars.at(i_plane))) {
	    ret = false;
	  }
	}
      }
      return ret;
    }
  }
}
