#include "RooStats/RooStatsUtils.h"
#include <iostream>
 
using namespace RooFit;
using namespace RooStats; // the utilities are in the RooStats namespace
using std::cout, std::endl;
 
void rs_nue_significance(const char *filename)
{
  // filename must be in this vector
  std::vector<std::string> filenames{"checkDataCuts_BDTcutFalse.root", "checkDataCuts_BDTcutTrue.root"};

  if (std::find(filenames.begin(), filenames.end(), filename) != filenames.end()) {
    // Filename is in the vector
    std::cout << "Valid filename: " << filename << std::endl;
  } else {
    // Filename is not in the vector
    std::cout << "Invalid filename: " << filename << std::endl;
    return
  }

  // WARNING WARNING WARNING SET TO 1 FOR REAL DATA ANALYSIS
  // Change only for sensitivity studies
  //float LUMI_SCALE = 2.5;
  float LUMI_SCALE = 1.;
  
  std::vector<std::string> mc_samples{"hadMC", "numuCC_MC", "nueCC_MC", "nuTauCC0mu_MC", "nuTauCC1mu_MC", "NC_MC"};
  // std::vector<std::string> mc_samples{"numuCC_MC", "nueCC_MC", "nuTauCC0mu_MC", "nuTauCC1mu_MC", "NC_MC"};

  std::map<std::string, std::string> sample_colors;
  sample_colors["hadMC"] = "#EE6677";
  sample_colors["numuCC_MC"] = "#22833";
  sample_colors["nueCC_MC"] = "#AA3377";
  sample_colors["NC_MC"] = "#CCBB44";
  sample_colors["nuTauCC0mu_MC"] = "#66CCEE";
  sample_colors["nuTauCC1mu_MC"] = "#4477AA";
  
  std::map<std::string, double> sample_scale;
  std::map<std::string, double> sample_fractional_uncertainty;

  // External constraints
  // --------------------------------------------------------------------------------
  //  sample_scale["numuCC_MC"] = (32-0.25)/19.1; // Based on extended numu analysis
  sample_scale["numuCC_MC"] = 1.; // Based on extended numu analysis
  sample_fractional_uncertainty["numuCC_MC"] = sqrt(32./pow(32,2) + pow(0.043,2)); // Statistical uncertainty on numu measurement added in quadrature with DS tracking systematic uncertainty (from flux measurement).

  sample_scale["nueCC_MC"] = 1.; // PLACEHOLDER
  sample_fractional_uncertainty["nueCC_MC"] = 1.; // PLACEHOLDER

  //  sample_scale["NC_MC"] = (32-0.25)/19.1; // Assume same scaling and uncertainty for NC. This is only used for the nueCC search.
  sample_scale["NC_MC"] = 1; // Assume same scaling and uncertainty for NC. This is only used for the nueCC search.
  sample_fractional_uncertainty["NC_MC"] = sqrt(32./pow(32,2) + pow(0.043,2)); //

  sample_scale["nuTauCC0mu_MC"] = 1.; // No external constraint
  sample_fractional_uncertainty["nuTauCC0mu_MC"] = 1.; // 100% uncertainty

  sample_scale["nuTauCC1mu_MC"] = 1.; // No external constraint
  sample_fractional_uncertainty["nuTauCC1mu_MC"] = 1.; // 100% uncertainty
  // --------------------------------------------------------------------------------
  
  auto f_histo = new TFile(filename);

  std::map<std::string, TH1D*> hists;
  for (std::string sample_name : mc_samples) hists[sample_name] = static_cast<TH1D*>(f_histo->Get(("hit_density_sel_"+sample_name).c_str()));
  for (std::string sample_name : mc_samples) hists[sample_name]->Scale(LUMI_SCALE);

  auto h_data = static_cast<TH1D*>(f_histo->Get("hit_density_sel"));
  h_data->Scale(LUMI_SCALE);
  
  // First, use background-dominated region to find the neutral hadron background scale and uncertainty:
  std::vector<double> had_CR_range{2000, 5000};
  std::vector<int> had_CR_bins{hists["hadMC"]->GetXaxis()->FindBin(had_CR_range.at(0)),hists["hadMC"]->GetXaxis()->FindBin(had_CR_range.at(1)-1)};

  std::cout << "CONTROL REGION EXPECTED BREAKDOWN" << std::endl;
  for (std::string sample_name : mc_samples){
    std::cout << sample_name << " " << hists[sample_name]->Integral(had_CR_bins.at(0), had_CR_bins.at(1)) << std::endl;
  }

  std::cout << "CONTROL REGION DATA" << std::endl;
  std::cout << h_data->Integral(had_CR_bins.at(0), had_CR_bins.at(1)) << std::endl;

  std::cout << "Neutral hadron background scale factor" << std::endl;
  double non_had_CR = 0.;
  for (std::string sample_name : mc_samples){
    if (sample_name == "hadMC") continue;

    non_had_CR += hists[sample_name]->Integral(had_CR_bins.at(0), had_CR_bins.at(1));
  }

  double had_CR = hists["hadMC"]->Integral(had_CR_bins.at(0), had_CR_bins.at(1));

  double data_CR = h_data->Integral(had_CR_bins.at(0), had_CR_bins.at(1));
  
  std::cout << "Total observed events in CR: " << data_CR << std::endl;
  std::cout << "Non-hadron events predicted in CR: " << non_had_CR << std::endl;
  std::cout << "Hadron events predicted in CR (factorised approach, large overestimation): " << had_CR << std::endl;

  double had_scale = (data_CR - non_had_CR)/had_CR;
  
  sample_scale["hadMC"] = had_scale;
  std::cout << "Hadron MC scale factor: " << sample_scale["hadMC"] << std::endl;
  sample_fractional_uncertainty["hadMC"] = 1.; // Set 100% uncertainty on the extrapolation

  // Now, find optimal region to exclude no-shower hypothesis using the MC expectation for the signal and the scaled background.
  auto g_SR_opt = new TGraph();
  auto g_SR_nueCC_opt = new TGraph();

  std::map<std::string, TGraphErrors*> g_SR_rates;

  for (std::string sample_name : mc_samples){
    g_SR_rates[sample_name] = new TGraphErrors();
  }

  for (double test_SR = 5000; test_SR <= 20000; test_SR += 500){
    std::cout << test_SR << std::endl;

    int bin = hists["hadMC"]->GetXaxis()->FindBin(test_SR);

    std::cout << bin << " " << test_SR << std::endl;
    
    for (std::string sample_name : mc_samples){

      double integral = hists[sample_name]->Integral(bin, hists[sample_name]->GetNbinsX()+1);

      g_SR_rates[sample_name]->AddPoint(test_SR, integral);

      if (sample_fractional_uncertainty.find(sample_name) != sample_fractional_uncertainty.end()) {
	g_SR_rates[sample_name]->SetPointError(g_SR_rates[sample_name]->GetN()-1, 0., sample_fractional_uncertainty[sample_name]);
      } else {
	g_SR_rates[sample_name]->SetPointError(g_SR_rates[sample_name]->GetN()-1, 0., 0.);
      }
    }
  }

  g_SR_rates["0mu"] = new TGraphErrors();

  for (int i_point = 0; i_point < g_SR_rates["NC_MC"]->GetN(); i_point++){

    double test_SR = g_SR_rates["NC_MC"]->GetPointX(i_point);
    
    g_SR_rates["0mu"]->AddPoint(test_SR,
				g_SR_rates["NC_MC"]->GetPointY(i_point) + g_SR_rates["nueCC_MC"]->GetPointY(i_point) + g_SR_rates["nuTauCC0mu_MC"]->GetPointY(i_point));

    double n_sig = 0.;
    double n_background = 0.;
    double uncert_background = 0.;

    double n_sig_nueCC = 0.;
    double n_background_nueCC = 0.;
    double uncert_background_nueCC = 0.;
    
    for (std::string sample_name : mc_samples){
      double rate = g_SR_rates[sample_name]->GetPointY(i_point);
      double error = g_SR_rates[sample_name]->GetErrorY(i_point);
      
      if (sample_name.compare("nueCC_MC") == 0) { // Signal in both cases
	n_sig_nueCC += rate*sample_scale[sample_name];
	n_sig += rate*sample_scale[sample_name];
      }
      else if ((sample_name.compare("NC_MC") == 0) or (sample_name.compare("nuTauCC0mu_MC") == 0)){ // Signal for 0mu, else background
	n_sig += rate*sample_scale[sample_name];
	
	n_background_nueCC += rate*sample_scale[sample_name];
	uncert_background_nueCC = sqrt(pow(uncert_background_nueCC, 2) + pow(error*sample_scale[sample_name]*rate, 2));
      } else { // Background in all cases

	n_background += rate*sample_scale[sample_name];
	uncert_background = sqrt(pow(uncert_background, 2) + pow(error*sample_scale[sample_name]*rate, 2));
	
	n_background_nueCC += rate*sample_scale[sample_name];
	uncert_background_nueCC = sqrt(pow(uncert_background_nueCC, 2) + pow(error*sample_scale[sample_name]*rate, 2));
      }
    }
    
    double significance = NumberCountingUtils::BinomialExpZ(n_sig, n_background, uncert_background/n_background);
    double significance_nueCC = NumberCountingUtils::BinomialExpZ(n_sig_nueCC, n_background_nueCC, uncert_background_nueCC/n_background_nueCC);
    
    std::cout <<" SR " << test_SR << " n_sig " << n_sig << " n_background " << n_background << " uncert_background " << uncert_background/n_background << " expected significance " << significance << std::endl; 
    std::cout << "NUE SR " << test_SR << " n_sig " << n_sig_nueCC << " n_background " << n_background_nueCC << " uncert_background " << uncert_background_nueCC/n_background_nueCC << " expected significance " << significance_nueCC << std::endl; 
    g_SR_opt->AddPoint(test_SR, significance);
    g_SR_nueCC_opt->AddPoint(test_SR, significance_nueCC);
  }

  auto c_opt = new TCanvas("c_opt");
  g_SR_opt->Draw("AL");
  g_SR_opt->GetXaxis()->SetTitle("Signal region SciFi hit density lower bound");
  g_SR_opt->GetYaxis()->SetTitle("Expected observation significance / #sigma");
  
  if ( strstr(filename, "True") != NULL ) {
    g_SR_opt->SetTitle("Expected observation significance with BDT cut");
  } else if ( strstr(filename, "False") != NULL ) {
    g_SR_opt->SetTitle("Expected observation significance without BDT cut");
  }

  g_SR_opt->GetYaxis()->SetRangeUser(0, 10.);
  
  g_SR_nueCC_opt->SetLineColor(kRed);
  g_SR_nueCC_opt->Draw("L");

  auto leg = new TLegend(0.6, 0.8, 0.9, 0.9);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(g_SR_opt, "#nu0#mu", "l");
  leg->AddEntry(g_SR_nueCC_opt, "#nu_{e}CC", "l");
  leg->Draw();
  gPad->Update();


  auto c_rates_0mu = new TCanvas("c_rates_0mu");
  auto leg_rates_0mu = new TLegend(0.6, 0.7, 0.9, 0.9);
  
  g_SR_rates["0mu"]->SetLineColor(TColor::GetColor(sample_colors["nueCC_MC"].c_str()));
  g_SR_rates["0mu"]->Draw("AL");
  g_SR_rates["0mu"]->GetYaxis()->SetRangeUser(0.0001, 80);
  g_SR_rates["0mu"]->GetYaxis()->SetTitle("Expected number of events in signal region");
  g_SR_rates["0mu"]->GetXaxis()->SetTitle("Signal region SciFi hit density lower bound");
  leg_rates_0mu->AddEntry(g_SR_rates["0mu"], "#nu0#mu", "lp");
  
  for (std::string sample_name : mc_samples){
    g_SR_rates[sample_name]->SetLineColor(TColor::GetColor(sample_colors[sample_name].c_str()));
    
    if ((sample_name.compare("NC_MC") == 0) or (sample_name.compare("nueCC_MC") == 0) or (sample_name.compare("nuTauCC0mu_MC") == 0)) continue;

    g_SR_rates[sample_name]->Draw("L");
    leg_rates_0mu->AddEntry(g_SR_rates[sample_name], sample_name.c_str(), "lp");
  }

  leg_rates_0mu->SetBorderSize(0);
  leg_rates_0mu->SetFillStyle(0);
  leg_rates_0mu->Draw();
  c_rates_0mu->SetLogy();
  c_rates_0mu->RedrawAxis();
  
  auto c_rates_nueCC = new TCanvas("c_rates_nueCC");
  auto leg_rates_nueCC = new TLegend(0.6, 0.7, 0.9, 0.9);
  bool first = true;
  for (std::string sample_name : mc_samples){
    g_SR_rates[sample_name]->SetLineColor(TColor::GetColor(sample_colors[sample_name].c_str()));
    leg_rates_nueCC->AddEntry(g_SR_rates[sample_name], sample_name.c_str(), "lp");
    if (first){
      g_SR_rates[sample_name]->Draw("AL");
      g_SR_rates[sample_name]->GetYaxis()->SetRangeUser(0.0001, 80);
      g_SR_rates[sample_name]->GetYaxis()->SetTitle("Expected number of events in signal region");
      g_SR_rates[sample_name]->GetXaxis()->SetTitle("Signal region SciFi hit density lower bound");
      first = false;
    } else {
      g_SR_rates[sample_name]->Draw("L");
    }
  }
  leg_rates_nueCC->SetBorderSize(0);
  leg_rates_nueCC->SetFillStyle(0);
  leg_rates_nueCC->Draw();
  c_rates_nueCC->SetLogy();
  c_rates_nueCC->RedrawAxis();

  // OBSERVED SIGNIFICANCE (all nutau background)
  double n_data = 9;

  // NOTE THAT THE NUMBERS BELOW ARE ALL HARD-CODED. COPY NUMBERS BY HAND TO OBTAIN OBSERVATION SIGNIFICANCE.
  // OBSERVED SIGNIFICANCE (nutau0mu signal)
  std::cout << "OBSERVED SIGNIFICANCE " << NumberCountingUtils::BinomialExpZ(n_data, 0.19518, 0.183765) << std::endl;
  // SR 11000 n_sig 6.35639 n_background 0.19518 uncert_background 0.183765 expected significance 5.49917

  std::cout << "Nue significance " << NumberCountingUtils::BinomialExpZ(n_data, 1.59674, 0.159669) << std::endl;
  // NUE SR 11000 n_sig 4.95624 n_background 1.59533 uncert_background 0.159755 expected significance 2.67008
  
  std::cout << "SIGNAL REGION EXPECTED BREAKDOWN" << std::endl;
  for (std::string sample_name : mc_samples){
    std::cout << sample_name << " " << hists[sample_name]->Integral(56 ,-1)*sample_scale[sample_name] << " +- " << sample_fractional_uncertainty[sample_name]*hists[sample_name]->Integral(56 ,-1)*sample_scale[sample_name] << std::endl;
  }
  std::cout << "ALL REGIONS EXPECTED BREAKDOWN" << std::endl;
  for (std::string sample_name : mc_samples){
    std::cout << sample_name << " " << hists[sample_name]->Integral(0 ,-1)*sample_scale[sample_name] << " +- " << sample_fractional_uncertainty[sample_name]*hists[sample_name]->Integral(56 ,-1)*sample_scale[sample_name] << std::endl;
  }

  
}
