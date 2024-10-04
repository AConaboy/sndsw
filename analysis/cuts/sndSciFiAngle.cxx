#include "sndSciFiAngle.h"

#include "sndSciFiTools.h"

#include "TChain.h"

#include "TVector3.h"

#include "TFitResult.h"

namespace snd{
  namespace analysis_cuts {
    scifiAngle::scifiAngle(double intercept, double slope, double max_chi2, TChain * ch) : sciFiBaseCut(ch){

      intercept_ = intercept;
      slope_ = slope;
      max_chi2_ = max_chi2;
      
      cutName = "SciFi angle";

      shortName = "SciFiAngle";
      nbins = std::vector<int>{100, 150, 100, 150};
      range_start = std::vector<double>{-2.5, 0, -2.5, 0};
      range_end = std::vector<double>{2.5, 150, 2.5, 150};
      plot_var = std::vector<double>{-1, -1, -1, -1};

      gv_ = new TGraph();
      gh_ = new TGraph();
    }
    
    bool scifiAngle::passCut(){
      initializeEvent();

      int n_planes_v = 5 - std::count(hits_per_plane_vertical.begin(), hits_per_plane_vertical.end(), 0);
      int n_planes_h = 5 - std::count(hits_per_plane_horizontal.begin(), hits_per_plane_horizontal.end(), 0);

      if (n_planes_v < 2 or n_planes_h < 2) {
	plot_var[0] = 1000;
	plot_var[1] = -1;
	plot_var[2] = 1000;
	plot_var[3] = -2;
	return false;
      }
      
      gv_->Set(0);
      gh_->Set(0);
      
      sndScifiHit * hit;
      TIter hitIterator(scifiDigiHitCollection);

      while ( (hit = (sndScifiHit*) hitIterator.Next()) ){
	if (hit->isValid()){

	  scifiDet->GetSiPMPosition(hit->GetDetectorID(), a_, b_);

	  if (hit->isVertical()){
	    gv_->SetPoint(gv_->GetN(), (a_.Z() + b_.Z())/2., (a_.X() + b_.X())/2.);
	  } else {
	    gh_->SetPoint(gh_->GetN(), (a_.Z() + b_.Z())/2., (a_.Y() + b_.Y())/2.);
	  }
	}
      }

      //      if (gv_->GetN() == 0 or gh_->GetN() == 0){
      //      }

      auto fitv = gv_->Fit("pol1", "SQN");
      double tanv = fitv->GetParams()[1];
      double reducedchi2v = fitv->Chi2()/fitv->Ndf();
      
      auto fith = gh_->Fit("pol1", "SQN");
      double tanh = fith->GetParams()[1];
      double reducedchi2h = fith->Chi2()/fith->Ndf();

      plot_var[0] = tanv;
      plot_var[1] = reducedchi2v;

      plot_var[2] = tanh;
      plot_var[3] = reducedchi2h;

      if (reducedchi2v <= intercept_ - (2.5-std::abs(tanv))*slope_) return false;
      if (reducedchi2h <= intercept_ - (2.5-std::abs(tanh))*slope_) return false;
      if (reducedchi2v > max_chi2_) return false;
      if (reducedchi2h > max_chi2_) return false;
      
      return true;
    }
  }	     
}
