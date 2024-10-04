#include "sndSciFiBaseCut.h"

#include <vector>

#include "TClonesArray.h"
#include "TChain.h"
#include "sndScifiHit.h"
#include "Scifi.h"
#include "TROOT.h"
#include "TH1D.h"
#include "SNDLHCEventHeader.h"

namespace snd::analysis_cuts {

  TChain * sciFiBaseCut::tree = 0;
  SNDLHCEventHeader * sciFiBaseCut::header = 0;
  unsigned long int sciFiBaseCut::read_entry = -1;
  
  TClonesArray * sciFiBaseCut::scifiDigiHitCollection = 0;
  TClonesArray * sciFiBaseCut::scifiDigiHitCollection_raw = 0;

  std::vector<TH1D *> sciFiBaseCut::hHitTime;
  
  std::vector<int> sciFiBaseCut::hits_per_plane_vertical = std::vector<int>(5, 0);
  std::vector<int> sciFiBaseCut::hits_per_plane_horizontal = std::vector<int>(5, 0);

  Scifi * sciFiBaseCut::scifiDet = 0;
  
  sciFiBaseCut::sciFiBaseCut(TChain * ch, bool select_events, double min_clock_cycle, double max_clock_cycle){

    if (tree == 0){
      tree = ch;

      if ( strcmp(tree->GetName(), "rawConv") == 0 ){
	// Real data times in clock cycle units
	TDC2ns_ = 1E9/160.316E6;
      } else {
	// MC in ns
	TDC2ns_ = 1.;
      }
      
      select_events_ = select_events;
      min_clock_cycle_ = min_clock_cycle;
      max_clock_cycle_ = max_clock_cycle;
      
      if (not select_events_) {
	scifiDigiHitCollection = new TClonesArray("sndScifiHit", 3000);
	tree->SetBranchAddress("Digi_ScifiHits", &scifiDigiHitCollection);
      }
      else {
	scifiDigiHitCollection_raw = new TClonesArray("sndScifiHit", 3000);
	tree->SetBranchAddress("Digi_ScifiHits", &scifiDigiHitCollection_raw);
	
	scifiDigiHitCollection = new TClonesArray("sndScifiHit", 3000);
	
	for (int i_orientation = 0; i_orientation < 2; i_orientation++){
	  for (int i_station = 0; i_station < 5; i_station++){
	    hHitTime.push_back(new TH1D(("hHitTime_"+std::to_string(i_orientation)+"_"+std::to_string(i_station)).c_str(), ";SciFi hit time [clock cycles]", 200., 0., 100.));
	  }
	}
      }
      
      scifiDet = (Scifi*) gROOT->GetListOfGlobals()->FindObject("Scifi");
	
      if (header == 0){
	// Check if branch exists in the TTree
	if (tree->FindBranch("EventHeader") != NULL){
	  if (tree->FindBranch("EventHeader")->GetAddress() != NULL) {
	    header = (SNDLHCEventHeader *) *((char**)tree->FindBranch("EventHeader")->GetAddress());
	    tree->GetEntry(0);
	    // Invalid runId. Invalid object?
	    if (header->GetRunId() == 0) header = NULL;
	  }
	}

	// Check alternate branch name
	if (tree->FindBranch("EventHeader.") != NULL and header == NULL){
	  if (tree->FindBranch("EventHeader.")->GetAddress() != NULL) {
	    header = (SNDLHCEventHeader *) *((char**) tree->FindBranch("EventHeader.")->GetAddress());
	    tree->GetEntry(0);
	    // Invalid runId. Invalid object?
	    if (header->GetRunId() == 0) header = NULL;
	  }
	}

	// Header is still null, create object and set address
	if (header == NULL) {
	  header = new SNDLHCEventHeader();
	  tree->SetBranchAddress("EventHeader", &header);
	  tree->GetEntry(0);
	  if (header->GetEventTime() == -1) {
	    tree->SetBranchAddress("EventHeader.", &header);
	    tree->GetEntry(0);
	    if (header->GetEventTime() == -1) throw std::runtime_error("Invalid event header");
	  }
	}
      }
      }
    }


  void sciFiBaseCut::initializeEvent(){
    if (read_entry != tree->GetReadEntry()){
      read_entry = tree->GetReadEntry();
      
      // Initialize detector
      scifiDet->InitEvent(header);

      sndScifiHit * hit;
      
      TIter hitIterator(scifiDigiHitCollection);
      TIter hitIterator_raw(scifiDigiHitCollection_raw);

      if (select_events_){
	for (int i_orientation = 0; i_orientation < 2; i_orientation++){
	  for (int i_station = 0; i_station < 5; i_station++){
	    hHitTime.at(i_orientation*5+i_station)->Reset();
	  }
	}
	scifiDigiHitCollection->Clear();
	
	while (hit = (sndScifiHit*) hitIterator_raw.Next()){
	  if (hit->isValid()){
	    int this_orientation = 0;
	    if (hit->isVertical()) this_orientation = 1;
	    int this_station = hit->GetStation() - 1;
	    // Get corrected time!!!
	    double time = hit->GetTime()*TDC2ns_;
	    time = scifiDet->GetCorrectedTime(hit->GetDetectorID(), time, 0);
	    hHitTime.at(this_orientation*5+this_station)->Fill(time);
	  }
	}

	std::vector<double> peakTiming(2*5, 0);

	for (int i_orientation = 0; i_orientation < 2; i_orientation++){
	  for (int i_station = 0; i_station < 5; i_station++){
	    int i = i_orientation*5+i_station;
	    peakTiming.at(i) = (hHitTime.at(i)->GetMaximumBin()-0.5)*(hHitTime.at(i)->GetXaxis()->GetXmax()-hHitTime.at(i)->GetXaxis()->GetXmin())/hHitTime.at(i)->GetNbinsX() + hHitTime.at(i)->GetXaxis()->GetXmin();
	  }
	}
	
	hitIterator_raw.Reset();
	int i_hit = 0;
	while (hit = (sndScifiHit*) hitIterator_raw.Next()){
	  if (not hit->isValid()) continue;
	  double time = hit->GetTime()*TDC2ns_;
	  time = scifiDet->GetCorrectedTime(hit->GetDetectorID(), time, 0);
	  int this_orientation = 0;
	  if (hit->isVertical()) this_orientation = 1;
	  int this_station = hit->GetStation() - 1;
	  int i = this_orientation*5+this_station;
	  // Cuts defined in clock cycle unit
	  if (time < peakTiming.at(i) + min_clock_cycle_*1E9/160.316E6) continue;
	  if (time > peakTiming.at(i) + max_clock_cycle_*1E9/160.316E6) continue;
	  (*scifiDigiHitCollection)[i_hit++] = hit;
	}
      }

      // Clear hits per plane vectors
      std::fill(hits_per_plane_vertical.begin(), hits_per_plane_vertical.end(), 0);
      std::fill(hits_per_plane_horizontal.begin(), hits_per_plane_horizontal.end(), 0);
      
      // Add valid hits to hits per plane vectors
      hitIterator.Reset();
      while ( (hit = (sndScifiHit*) hitIterator.Next()) ){
	if (hit->isValid()){
	  int sta = hit->GetStation();
	  if (hit->isVertical()){
	    hits_per_plane_vertical[sta-1]++;
	  } else {
	    hits_per_plane_horizontal[sta-1]++;
	  }
	}
      }
    }
  }
  sciFiBaseCut::~sciFiBaseCut(){
    delete scifiDigiHitCollection;
    delete scifiDigiHitCollection_raw;
    header = NULL;
  };
}
