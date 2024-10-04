#include <iostream>
#include <vector>

#include "TObject.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "TH1D.h"
#include "TPython.h"

#include "ShipMCTrack.h"

// Cuts
#include "sndBaseCut.h"
#include "sndMinSciFiHitsCut.h"
#include "sndSciFiStationCut.h"
#include "sndVetoCut.h"
#include "sndMinSciFiConsecutivePlanes.h"
#include "sndDSActivityCut.h"
#include "sndUSQDCCut.h"
#include "sndEventDeltat.h"
#include "sndAvgSciFiFiducialCut.h"
#include "sndAvgDSFiducialCut.h"
#include "sndDSVetoCut.h"
#include "sndSciFiContinuity.h"
#include "sndUSPlanesHit.h"
#include "sndUSBarsVeto.h"
#include "sndSciFiAngle.h"
#include "sndStableBeams.h"
#include "sndIP1.h"
#include "sndPrescale.h"

// Alternatice sets of cuts.
enum Cutset { stage1cuts, novetocuts, FVsideband, allowWalls2and5, stage1cutsVetoFirst, nueFilter, allowWalls2and5VetoSideband} ;

int main(int argc, char ** argv) {

  std::cout << "Starting neutrino filter" << std::endl;
  
  if (argc != 4) {
    std::cout << "Two arguments required: input file name (or reg exp), output file name, cut set (0: stage 1 selection, 1: no veto or scifi 1st layer, 2: FV sideband, 3: include walls 2 and 5, 4: apply veto cut first, 5: nue filter (no DS selection, allow walls 1 to 4." << std::endl;
    exit(-1);
  }

  // SNDSW geometry:
  TPython::Exec("import SndlhcGeo");
  // REMOVE REMOVE REMOVE hard-coded path!!!
  TPython::Exec("sndgeo = SndlhcGeo.GeoInterface('/eos/experiment/sndlhc/convertedData/physics/2023/geofile_sndlhc_TI18_V4_2023.root')");
  // Hacky hacky
  TPython::Exec("import atexit");
  TPython::Exec("def pyExit():\n print(\"Make suicide until solution found for freezing\")\n os.system('kill '+str(os.getpid()))\natexit.register(pyExit)");
  
  // Input files
  bool isMC = false;
  TChain * ch = new TChain("rawConv");
  ch->Add(argv[1]);
  if (ch->GetEntries() == 0){
    delete ch;
    ch = new TChain("cbmsim");
    ch->Add(argv[1]);
    if (ch->GetEntries() > 0) {
      isMC = true;
      ch->SetBranchStatus("EventHeader", 0); // To be able to run on old MC.
      ch->SetBranchStatus("EventHeader.", 0); // To be able to run on old MC. Bizarre....?
    } else {
      std::cout << "Didn't find rawConv or cbmsim in input file" << std::endl;
      exit(-1);
    }
  }
  std::cout << "Got input tree" << std::endl;

  // MC truth
  TClonesArray * MCTracks = new TClonesArray("ShipMCTrack", 5000);
  if (isMC) ch->SetBranchAddress("MCTrack", &MCTracks);

  // Output file
  TFile * outFile = new TFile(argv[2], "RECREATE");
  std::cout << "Got output file" << std::endl;  

  ch->GetEntry(0);
  ch->GetFile()->Get("BranchList")->Write("BranchList", TObject::kSingleKey);
  ch->GetFile()->Get("TimeBasedBranchList")->Write("TimeBasedBranchList", TObject::kSingleKey);
  if (ch->GetFile()->Get("FileHeader")) ch->GetFile()->Get("FileHeader")->Write();
  if (ch->GetFile()->Get("FileHeaderHeader")) ch->GetFile()->Get("FileHeaderHeader")->Write();

  // Set up all branches to copy to output TTree.
  TTree * outTree = ch->CloneTree(0);
  std::cout << "Got output tree" << std::endl;  

  // Set up cuts
  std::cout << "Starting cut set up" << std::endl;

  std::vector< snd::analysis_cuts::baseCut * > cutFlow;
  
  int selected_cutset = std::atoi(argv[3]);

  if (selected_cutset == stage1cuts){ // Stage 1 cuts
    if (not isMC) {
      cutFlow.push_back( new snd::analysis_cuts::stableBeamsCut(ch));
      cutFlow.push_back( new snd::analysis_cuts::ip1Cut(ch));
      cutFlow.push_back( new snd::analysis_cuts::eventDeltatCut(-1, 100, ch)); // J. Previous event more than 100 clock cycles away. To avoid deadtime issues.
    }
    cutFlow.push_back( new snd::analysis_cuts::avgSciFiFiducialCut(200, 1200, 300, 128*12-200, ch)); // E. Average SciFi hit channel number must be within [200, 1200] (ver) and [300, max-200] (hor)
    cutFlow.push_back( new snd::analysis_cuts::avgDSFiducialCut(70, 105, 10, 50, ch)); // F. Average DS hit bar number must be within [70, 105] (ver) and [10, 50] (hor)
    cutFlow.push_back( new snd::analysis_cuts::vetoCut(ch)); // B. No veto hits
    cutFlow.push_back( new snd::analysis_cuts::sciFiStationCut(0., std::vector<int>(1, 1), ch)); // C. No hits in first SciFi plane
    cutFlow.push_back( new snd::analysis_cuts::sciFiStationCut(0., std::vector<int>(1, 2), ch)); // C. No hits in second SciFi plane
    cutFlow.push_back( new snd::analysis_cuts::sciFiStationCut(0.05, std::vector<int>(1, 5), ch)); // D. Vertex not in 5th wall
    cutFlow.push_back( new snd::analysis_cuts::minSciFiConsecutivePlanes(ch)); // G. At least two consecutive SciFi planes hit
    cutFlow.push_back( new snd::analysis_cuts::DSActivityCut(ch)); // H. If there is a downstream hit, require hits in all upstream stations.    
  } else if (selected_cutset == stage1cutsVetoFirst){ // Stage 1 cuts but with Veto cut upfront. For neutral hadron background estimation
    if (not isMC) {
      cutFlow.push_back( new snd::analysis_cuts::stableBeamsCut(ch));
      cutFlow.push_back( new snd::analysis_cuts::ip1Cut(ch));
      cutFlow.push_back( new snd::analysis_cuts::eventDeltatCut(-1, 100, ch)); // J. Previous event more than 100 clock cycles away. To avoid deadtime issues.
    }
    cutFlow.push_back( new snd::analysis_cuts::vetoCut(ch)); // B. No veto hits
    cutFlow.push_back( new snd::analysis_cuts::avgSciFiFiducialCut(200, 1200, 300, 128*12-200, ch)); // E. Average SciFi hit channel number must be within [200, 1200] (ver) and [300, max-200] (hor)
    cutFlow.push_back( new snd::analysis_cuts::avgDSFiducialCut(70, 105, 10, 50, ch)); // F. Average DS hit bar number must be within [70, 105] (ver) and [10, 50] (hor)
    cutFlow.push_back( new snd::analysis_cuts::sciFiStationCut(0., std::vector<int>(1, 1), ch)); // C. No hits in first SciFi plane
    cutFlow.push_back( new snd::analysis_cuts::sciFiStationCut(0., std::vector<int>(1, 2), ch)); // C. No hits in second SciFi plane
    cutFlow.push_back( new snd::analysis_cuts::sciFiStationCut(0.05, std::vector<int>(1, 5), ch)); // D. Vertex not in 5th wall
    cutFlow.push_back( new snd::analysis_cuts::minSciFiConsecutivePlanes(ch)); // G. At least two consecutive SciFi planes hit
    cutFlow.push_back( new snd::analysis_cuts::DSActivityCut(ch)); // H. If there is a downstream hit, require hits in all upstream stations.    
  } else if (selected_cutset == novetocuts) {
    if (not isMC) {
      cutFlow.push_back( new snd::analysis_cuts::stableBeamsCut(ch));
      cutFlow.push_back( new snd::analysis_cuts::ip1Cut(ch));
      cutFlow.push_back( new snd::analysis_cuts::eventDeltatCut(-1, 100, ch)); // J. Previous event more than 100 clock cycles away. To avoid deadtime issues.
    }
    cutFlow.push_back( new snd::analysis_cuts::sciFiStationCut(0.05, std::vector<int>(1, 5), ch)); // D. Vertex not in 5th wall
    cutFlow.push_back( new snd::analysis_cuts::avgSciFiFiducialCut(200, 1200, 300, 128*12-200, ch)); // E. Average SciFi hit channel number must be within [200, 1200] (ver) and [300, max-200] (hor)
    cutFlow.push_back( new snd::analysis_cuts::avgDSFiducialCut(70, 105, 10, 50, ch)); // F. Average DS hit bar number must be within [70, 105] (ver) and [10, 50] (hor)
    cutFlow.push_back( new snd::analysis_cuts::minSciFiConsecutivePlanes(ch)); // G. At least two consecutive SciFi planes hit
    cutFlow.push_back( new snd::analysis_cuts::DSActivityCut(ch)); // H. If there is a downstream hit, require hits in all upstream stations.    
  } else if (selected_cutset == FVsideband){
    if (not isMC) {
      cutFlow.push_back( new snd::analysis_cuts::stableBeamsCut(ch));
      cutFlow.push_back( new snd::analysis_cuts::ip1Cut(ch));
      cutFlow.push_back( new snd::analysis_cuts::eventDeltatCut(-1, 100, ch)); // J. Previous event more than 100 clock cycles away. To avoid deadtime issues.
    }
    cutFlow.push_back( new snd::analysis_cuts::vetoCut(ch)); // B. No veto hits
    cutFlow.push_back( new snd::analysis_cuts::sciFiStationCut(0., std::vector<int>(1, 1), ch)); // C. No hits in first SciFi plane
    cutFlow.push_back( new snd::analysis_cuts::sciFiStationCut(0., std::vector<int>(1, 2), ch)); // D. Vertex not in 5th wall
    cutFlow.push_back( new snd::analysis_cuts::sciFiStationCut(0.05, std::vector<int>(1, 5), ch)); // D. Vertex not in 5th wall
    cutFlow.push_back( new snd::analysis_cuts::avgSciFiFiducialCut(200, 1200, 300, 128*12-200, ch, true)); // E. Average SciFi hit channel number must be within [200, 1200] (ver) and [300, max-200] (hor)
    cutFlow.push_back( new snd::analysis_cuts::minSciFiConsecutivePlanes(ch)); // G. At least two consecutive SciFi planes hit
    cutFlow.push_back( new snd::analysis_cuts::DSActivityCut(ch)); // H. If there is a downstream hit, require hits in all upstream stations.    
  } else if (selected_cutset == allowWalls2and5) {
    if (not isMC) {
      cutFlow.push_back( new snd::analysis_cuts::stableBeamsCut(ch));
      cutFlow.push_back( new snd::analysis_cuts::ip1Cut(ch));
      cutFlow.push_back( new snd::analysis_cuts::eventDeltatCut(-1, 100, ch)); // J. Previous event more than 100 clock cycles away. To avoid deadtime issues.
    }
    cutFlow.push_back( new snd::analysis_cuts::avgSciFiFiducialCut(200, 1200, 300, 128*12-200, ch)); // E. Average SciFi hit channel number must be within [200, 1200] (ver) and [300, max-200] (hor)
    cutFlow.push_back( new snd::analysis_cuts::avgDSFiducialCut(70, 105, 10, 50, ch)); // F. Average DS hit bar number must be within [70, 105] (ver) and [10, 50] (hor)
    cutFlow.push_back( new snd::analysis_cuts::vetoCut(ch)); // B. No veto hits
    cutFlow.push_back( new snd::analysis_cuts::sciFiStationCut(0., std::vector<int>(1, 1), ch)); // C. No hits in first SciFi plane
    cutFlow.push_back( new snd::analysis_cuts::DSActivityCut(ch)); // H. If there is a downstream hit, require hits in all upstream stations.    
  } else if (selected_cutset == allowWalls2and5VetoSideband) {
    if (not isMC) {
      cutFlow.push_back( new snd::analysis_cuts::preScale(0.001) ); // Save only 1 in 1000 events.
      cutFlow.push_back( new snd::analysis_cuts::stableBeamsCut(ch));
      cutFlow.push_back( new snd::analysis_cuts::ip1Cut(ch));
      cutFlow.push_back( new snd::analysis_cuts::eventDeltatCut(-1, 100, ch)); // J. Previous event more than 100 clock cycles away. To avoid deadtime issues.
    }
    cutFlow.push_back( new snd::analysis_cuts::avgSciFiFiducialCut(200, 1200, 300, 128*12-200, ch)); // E. Average SciFi hit channel number must be within [200, 1200] (ver) and [300, max-200] (hor)
    cutFlow.push_back( new snd::analysis_cuts::avgDSFiducialCut(70, 105, 10, 50, ch)); // F. Average DS hit bar number must be within [70, 105] (ver) and [10, 50] (hor)
    cutFlow.push_back( new snd::analysis_cuts::DSActivityCut(ch)); // H. If there is a downstream hit, require hits in all upstream stations.    
  } else if (selected_cutset == nueFilter) {
    if (not isMC) {
      cutFlow.push_back( new snd::analysis_cuts::stableBeamsCut(ch));
      cutFlow.push_back( new snd::analysis_cuts::ip1Cut(ch));
      cutFlow.push_back( new snd::analysis_cuts::eventDeltatCut(-1, 100, ch)); // J. Previous event more than 100 clock cycles away. To avoid deadtime issues.
    }
    cutFlow.push_back( new snd::analysis_cuts::avgSciFiFiducialCut(200, 1200, 300, 128*12-200, ch)); // E. Average SciFi hit channel number must be within [200, 1200] (ver) and [300, max-200] (hor)
    cutFlow.push_back( new snd::analysis_cuts::USBarsVeto(std::vector<std::pair<int, double> >{{0, 2.}, {1, 2.}}, true, ch)); // Reject events with hits in lowest bars in first two US planes
    cutFlow.push_back( new snd::analysis_cuts::USBarsVeto(std::vector<std::pair<int, double> >{{0, 8.}, {1, 8.}}, false, ch)); // Reject events with hits in highest bars in first two US planes
    cutFlow.push_back( new snd::analysis_cuts::vetoCut(ch)); // B. No veto hits
    cutFlow.push_back( new snd::analysis_cuts::minSciFiConsecutivePlanes(ch)); // G. At least two consecutive SciFi planes hit
    cutFlow.push_back( new snd::analysis_cuts::sciFiContinuity(ch)); // All SciFi planes downstream of first active (both views) plane must be hit (both views).
    cutFlow.push_back( new snd::analysis_cuts::USPlanesHit(std::vector<int>{0, 1}, ch)); // All SciFi planes downstream of first active (both views) plane must be hit (both views).
    cutFlow.push_back( new snd::analysis_cuts::minSciFiHits(35, ch)); // At least 35 SciFi hits
    if (isMC) cutFlow.push_back( new snd::analysis_cuts::USQDCCut(700, ch)); // Min QDC
    else      cutFlow.push_back( new snd::analysis_cuts::USQDCCut(600, ch)); // 
    cutFlow.push_back( new snd::analysis_cuts::DSVetoCut(ch)); // D. Veto events with hits in last DS planes
  } else {
    std::cout << "Unrecognized cutset. Exitting" << std::endl;
    exit(-1);
  }

  std::cout << "Done initializing cuts" << std::endl;

  int n_cuts = (int) cutFlow.size();

  // Book histograms
  // Cut-by-cut
  // All cut variables
  std::vector<std::vector<TH1D*> > cut_by_cut_var_histos = std::vector<std::vector<TH1D*> >();
  for (int i_cut = -1; i_cut < n_cuts; i_cut++){
    std::vector<TH1D*> this_cut_by_cut_var_histos = std::vector<TH1D*>();
    for (snd::analysis_cuts::baseCut * cut : cutFlow) {
      for(int i_dim = 0; i_dim < cut->getNbins().size(); i_dim++){
	this_cut_by_cut_var_histos.push_back(new TH1D((std::to_string(i_cut)+"_"+cut->getShortName()+"_"+std::to_string(i_dim)).c_str(), 
						      cut->getShortName().c_str(), 
						      cut->getNbins()[i_dim], cut->getRangeStart()[i_dim], cut->getRangeEnd()[i_dim]));
      }
    }
    cut_by_cut_var_histos.push_back(this_cut_by_cut_var_histos);
  }
  
  std::vector<std::vector<std::vector<TH1D*> > > cut_by_cut_truth_histos = std::vector<std::vector<std::vector<TH1D*> > >();
  if (isMC) {
    for (int i_species = 0; i_species < 6; i_species++){ // e, mu, tau0mu, tau1mu, NC, other
      std::vector<std::vector<TH1D*> > this_species_histos = std::vector<std::vector<TH1D*> >();
      std::string species_suffix;
      switch (i_species) {
      case 0:
	species_suffix = "nueCC";
	break;
      case 1:
	species_suffix = "numuCC";
	break;
      case 2:
	species_suffix = "nutauCC0mu";
	break;
      case 3:
	species_suffix = "nutauCC1mu";
	break;
      case 4:
	species_suffix = "NC";
	break;
      case 5:
	species_suffix = "Other"; // For PG sim
	break;
      default :
	std::cerr << "MC truth histograms initialization error! Unknown species" << std::endl;
	exit(-1);
      }

      for (int i_cut = -1; i_cut < n_cuts; i_cut++){
	
	std::vector<TH1D*> this_cut_by_cut_truth_histos = std::vector<TH1D*>();
	this_cut_by_cut_truth_histos.push_back(new TH1D((species_suffix+"_"+std::to_string(i_cut)+"_Enu").c_str(), "Enu", 300, 0, 3000));
	this_cut_by_cut_truth_histos.push_back(new TH1D((species_suffix+"_"+std::to_string(i_cut)+"_EEM").c_str(), "ELep", 300, 0, 3000));
	this_cut_by_cut_truth_histos.push_back(new TH1D((species_suffix+"_"+std::to_string(i_cut)+"_EHad").c_str(), "EHad", 300, 0, 3000));
	this_cut_by_cut_truth_histos.push_back(new TH1D((species_suffix+"_"+std::to_string(i_cut)+"_vtxX").c_str(), "vtxX", 200, -100, 0));
	this_cut_by_cut_truth_histos.push_back(new TH1D((species_suffix+"_"+std::to_string(i_cut)+"_vtxY").c_str(), "vtxY", 200, 0, 100));
	this_cut_by_cut_truth_histos.push_back(new TH1D((species_suffix+"_"+std::to_string(i_cut)+"_vtxZ").c_str(), "vtxZ", 200, 280, 380));

	this_species_histos.push_back(this_cut_by_cut_truth_histos);
      }
      cut_by_cut_truth_histos.push_back(this_species_histos);
    }
  }
  // N-1
  std::vector<TH1D*> n_minus_1_var_histos = std::vector<TH1D*>();
  for (snd::analysis_cuts::baseCut * cut : cutFlow) {
    for(int i_dim = 0; i_dim < cut->getNbins().size(); i_dim++){
      n_minus_1_var_histos.push_back(new TH1D(("n_minus_1_"+cut->getShortName()+"_"+std::to_string(i_dim)).c_str(), 
					      cut->getShortName().c_str(), 
					      cut->getNbins()[i_dim], cut->getRangeStart()[i_dim], cut->getRangeEnd()[i_dim]));
    }
  }

  std::cout << "Done initializing histograms" << std::endl;

  // Cut flow
  TH1D * cutFlowHistogram = new TH1D("cutFlow", "Cut flow;;Number of events passing cut", n_cuts+1, 0, n_cuts+1);
  for (int i = 2; i <= cutFlowHistogram->GetNbinsX(); i++){
    cutFlowHistogram->GetXaxis()->SetBinLabel(i, cutFlow.at(i-2)->getName().c_str());
  }

  std::cout << "Done initializing cut flow histogram" << std::endl;

  // Get number of entries
  unsigned long int n_entries = ch->GetEntries();

  // Holder for cut results
  std::vector<bool> passes_cut = std::vector<bool>(n_cuts, false);

  std::cout << "Starting event loop" << std::endl;
  for (unsigned long int i_entry = 0; i_entry < n_entries; i_entry++){
    ch->GetEntry(i_entry);
    if (i_entry%10000 == 0) std::cout << "Reading entry " << i_entry << std::endl;
    
    cutFlowHistogram->Fill(0);

    // Apply cuts
    int n_cuts_passed = 0;
    bool accept_event = true;
    int i_cut = 0;
    for (snd::analysis_cuts::baseCut * cut : cutFlow){
      if (cut->passCut()){
	if (accept_event) cutFlowHistogram->Fill(1 + i_cut);
	passes_cut[i_cut] = true;
	n_cuts_passed++;
      } else {
	accept_event = false;
	passes_cut[i_cut] = false;
      }
      i_cut++;
    }
    if (accept_event) outTree->Fill();


    // Fill histograms
    std::vector<TH1D*>::iterator hist_it;
    // Sequential
    for (int seq_cut = -1; seq_cut < ((int) passes_cut.size()); seq_cut++){
      if (seq_cut >= 0){
	if (not passes_cut[seq_cut]) break;
      }
      hist_it = cut_by_cut_var_histos[seq_cut+1].begin();
      for (snd::analysis_cuts::baseCut * cut : cutFlow) {
	for (int i_dim = 0; i_dim < cut->getPlotVar().size(); i_dim++){
	  (*hist_it)->Fill(cut->getPlotVar()[i_dim]);
	  hist_it++;
	}
      }
      if (isMC) {
	
	int this_species = -1;
	if (MCTracks->GetEntries() < 2) {
	  this_species = 5;
	} else {
	  
	  int pdgIn = abs(((ShipMCTrack*)MCTracks->At(0))->GetPdgCode());
	  int pdgOut = abs(((ShipMCTrack*)MCTracks->At(1))->GetPdgCode());
	
	  if (pdgIn == (pdgOut+1)){
	    //CC
	    if (pdgIn == 12) this_species = 0; // nueCC
	    if (pdgIn == 14) this_species = 1; // numuCC
	    if (pdgIn == 16) {
	      bool is1Mu = false;
	      for(int j_track = 2; j_track < MCTracks->GetEntries(); j_track++){
		if (((ShipMCTrack*)MCTracks->At(j_track))->GetMotherId() == 1 and abs(((ShipMCTrack*)MCTracks->At(j_track))->GetPdgCode()) == 13){
		  is1Mu = true;
		  break;
		}
	      }
	      if (is1Mu) this_species = 3; // nutauCC1mu
	      else this_species = 2; // nutauCC0mu
	    }
	  } else if (pdgIn == pdgOut) {
	    //NC
	    this_species = 4;
	  } else {
	    // Other
	    this_species = 5;
	  }
	}
	cut_by_cut_truth_histos[this_species][seq_cut+1][0]->Fill(((ShipMCTrack*) MCTracks->At(0))->GetEnergy()); // Enu
	if (this_species < 5) {
	  cut_by_cut_truth_histos[this_species][seq_cut+1][1]->Fill(((ShipMCTrack*) MCTracks->At(1))->GetEnergy()); // ELep
	  cut_by_cut_truth_histos[this_species][seq_cut+1][2]->Fill(((ShipMCTrack*) MCTracks->At(0))->GetEnergy()-((ShipMCTrack*) MCTracks->At(1))->GetEnergy()); // EHad
	}
	cut_by_cut_truth_histos[this_species][seq_cut+1][3]->Fill(((ShipMCTrack*) MCTracks->At(0))->GetStartX()); // X
	cut_by_cut_truth_histos[this_species][seq_cut+1][4]->Fill(((ShipMCTrack*) MCTracks->At(0))->GetStartY()); // Y
	cut_by_cut_truth_histos[this_species][seq_cut+1][5]->Fill(((ShipMCTrack*) MCTracks->At(0))->GetStartZ()); // Z
      }
    }

    // N-1
    int current_cut = 0;
    hist_it = n_minus_1_var_histos.begin();
    for (snd::analysis_cuts::baseCut * cut : cutFlow) {
      for (int i_dim = 0; i_dim < cut->getPlotVar().size(); i_dim++){
	if (((not passes_cut[current_cut]) and (n_cuts_passed == (cutFlow.size()-1)))
	    or (n_cuts_passed == cutFlow.size()) ) (*hist_it)->Fill(cut->getPlotVar()[i_dim]);
	hist_it++;
      }
      current_cut++;
    }
  }
  
  outFile->Write();

  delete outTree;
  delete outFile;
  delete MCTracks;
  for (auto cut : cutFlow) delete cut;
  delete ch;
  
  TPython::Exec("del sndgeo");
  
  return 0;
}
