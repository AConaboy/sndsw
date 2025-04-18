import ROOT
import os
from array import array
import random
from pathlib import Path
from tqdm import tqdm
from argparse import ArgumentParser

LUMI = 68.551

# SCALE TOTAL HADRON CONTRIBUTION TO THIS NUMBER. THIS IS the 90% upper limit OBTAINED FROM muonDISlumiCalc.py)
N_HAD = 8.21

BASE_FILTERED_DIR = Path("/eos/user/c/cvilela/SND_nue_analysis_May24/")

ch = ROOT.TChain("rawConv")

import SndlhcGeo
snd_geo = SndlhcGeo.GeoInterface("/afs/cern.ch/user/a/aconsnd/geofile_full.Ntuple-TGeant4-0.root") # geofile that contains the most up to date sim parameters
scifiDet = ROOT.gROOT.GetListOfGlobals().FindObject('Scifi')
muFilterDet = ROOT.gROOT.GetListOfGlobals().FindObject('MuFilter')

from sciFiTools import *

"""
Hacky way of determining BDT features until barycentre determination is incorportated into the mufilter & mufilter hit code

1. Make sure simulation flag changes
2. Load in correction and alignment parameters for the data
3. Make sure the alignment parameters change depending on the time alignment
"""

import Monitor, TimeWalk, SndlhcTracking, joblib, json
from args_config import add_arguments
from AnalysisFunctions import Analysis

parser=ArgumentParser()
parser.add_argument("--BDTcut", dest="BDTcut", help="bool for applying BDT cut", action='store_true')
add_arguments(parser)
options = parser.parse_args()

options.path='/eos/experiment/sndlhc/convertedData/physics/2023/'
options.mode='extendedreconstruction'
options.referencesystem=1 # Use Scifi alignment

muAna = Analysis(options)

from HCALTools import HCALTools
import xgboost as xgb
hcalTools = HCALTools(muAna, muFilterDet)
hcalTools.filekey=0

# Load in trained BDT
hcalTools.model = xgb.XGBClassifier(use_label_encoder=False, eval_metric='logloss', random_state=42)
hcalTools.model.load_model('/eos/home-a/aconsnd/SWAN_projects/Data analysis/bdt_model.json')

from numpy import load
hcalTools.BDTthreshold = load('/eos/home-a/aconsnd/SWAN_projects/Data analysis/best_threshold.npy')

# Load in the features used to train the bdt
# with open('/eos/home-a/aconsnd/SWAN_projects/Data analysis/bdt_features.json') as jf:
#     hcalTools.BDT_features = json.load(jf)

# Tell muAna to look for reference time in hcalTools
muAna.SetTask(hcalTools)

for this_run in (BASE_FILTERED_DIR / "data_2022_2023").glob("*/filtered_*.root"):
    ch.Add(this_run.as_posix())

N_MC_FILES=400

"""
For the Monte Carlo, use the MC that has been redigitised with the new timing digitisation that is compatible with data.
"""

chMC = ROOT.TChain("cbmsim")
# chMC.Add((BASE_FILTERED_DIR / "nuMC" / "filtered_stage1.root").as_posix())
redigitisedMC_path = '/eos/experiment/sndlhc/users/aconsnd/simulation/neutrino/data/sndlhc_13TeV_down_volTarget_100fb-1_SNDG18_02a_01_000/nueFilter/'
allfiles=os.listdir(redigitisedMC_path)
for fname in allfiles:
    if fname.endswith('root'): chMC.Add(redigitisedMC_path+fname)

chNeutral = ROOT.TChain("cbmsim")
chNeutral.Add("/afs/cern.ch/work/c/cvilela/public/SND_Nov_2023/sndsw/analysis/scripts/neutron_kaon_nue_stage1_noprescale.root")

def makePlots(ch, BDT_cut, name = "", isNuMC = False, 
            isNeutralHad = False, preselection = 1.):
    n = [0]*5

    selected_list = []
    tdiff_list = []
    max_slope_list = []
    
    h_n_hits = []
    h_n_hits_sel = []
    h_hit_density = []
    h_hit_density_sel = []
    h_hit_density_sel_after_dens2 = []
    h_hit_density_sel_after_dens = []
    h_hit_density_sel_precut = []
    h_log_hit_density_sel = []
    h_hit_density2_sel_precut = []
    h_hit_density2_sel = []
    
    h_min_chi2 = []
    h_log_min_chi2 = []
    h_SciFiAngle = []
    h_SciFiAngle_v_chi2 = []
    h_SciFiAngle_h_chi2 = []
    
    h_theta = []
    h_theta_density = []

    if isNuMC:
        flav_list = ["_NC", "_nueCC", "_numuCC", "_nuTauCC0mu", "_nuTauCC1mu"]
    else:
        flav_list = [""]

    BDT_weights = []

    for flav in flav_list:

        # Histograms for BDT signal probabilities
        BDT_weights.append(ROOT.TH1D("BDT_signalProb"+flav+name, "Prob event has no muon", 100, 0.0, 1.0 + 1e-6))

        h_n_hits.append(ROOT.TH1D("n_hits"+flav+name, ";Number of in-time SciFi hits", 1000, 0, 5000))
        h_n_hits_sel.append(ROOT.TH1D("n_hits_sel"+flav+name, ";Number of in-time SciFi hits", 1000, 0, 5000))
        h_hit_density.append(ROOT.TH1D("hit_density"+flav+name, ";Sum of hit densities", 200, 0, 40000))
        h_hit_density_sel.append(ROOT.TH1D("hit_density_sel"+flav+name, ";Sum of hit densities", 200, 0, 40000))
        h_hit_density_sel_precut.append(ROOT.TH1D("hit_density_sel_precut"+flav+name, ";Sum of hit densities", 200, 0, 40000))
        h_hit_density_sel_after_dens2.append(ROOT.TH1D("hit_density_sel_after_dens2"+flav+name, ";Sum of hit densities", 200, 0, 40000))
        h_hit_density_sel_after_dens.append(ROOT.TH1D("hit_density_sel_after_dens"+flav+name, ";Sum of hit densities", 200, 0, 40000))
        h_log_hit_density_sel.append(ROOT.TH1D("log_hit_density_sel"+flav+name, ";log(Sum of hit densities)", 120, 0, 12))

        h_hit_density2_sel_precut.append(ROOT.TH1D("hit_density2_sel_precut"+flav+name, ";Sum of hit densities", 500, 0, 10000))
        h_hit_density2_sel.append(ROOT.TH1D("hit_density2_sel"+flav+name, ";Sum of hit densities", 500, 0, 10000))

        h_min_chi2.append(ROOT.TH1D("min_chi2"+flav+name, ";min chi2", 100, 0., 10000.))
        h_log_min_chi2.append(ROOT.TH1D("log_min_chi2"+flav+name, ";log(min chi2)", 120, 0., 12.))
        h_SciFiAngle.append(ROOT.TH2D("SciFiAngle"+flav+name, ";SciFiAngle v;SciFiAngle h", 100, -2.5, 2.5, 100, -2.5, 2.5))
        h_SciFiAngle_v_chi2.append(ROOT.TH2D("SciFiAngle_v_chi2"+flav+name, ";SciFiAngle v;chi2", 100, -2.5, 2.5, 150, 0., 150.))
        h_SciFiAngle_h_chi2.append(ROOT.TH2D("SciFiAngle_h_chi2"+flav+name, ";SciFiAngle h;chi2", 100, -2.5, 2.5, 150, 0., 150.))
            
        h_theta.append(ROOT.TH1D("theta"+flav+name, ";#theta", 100, 0, 1))
        h_theta_density.append(ROOT.TH2D("theta_density"+flav+name, ";#theta;Sum of hit densities", 100, 0, 1, 200, 0, 40000))

    BDT_weights.append(ROOT.TH1D("BDT_signalProb_all", "Prob event has no muon", 100, 0.0, 1.0 + 1e-6))

    if isNeutralHad:
        totWeight = 0
    
    for i_event, event in tqdm(enumerate(ch)):

        # print(f'event {i_event}')

        if preselection < 1.:
            if random.random() > preselection:
                continue
            
        if not (isNuMC or isNeutralHad):
            scifiDet.InitEvent(event.EventHeader)
            muFilterDet.InitEvent(event.EventHeader)

        i_flav = 0
        if isNuMC:
            if event.MCTrack[0].GetPdgCode() == event.MCTrack[1].GetPdgCode():
                i_flav = 0 #NC
            elif abs(event.MCTrack[1].GetPdgCode()) == 11:
                i_flav = 1 #nueCC
            elif abs(event.MCTrack[1].GetPdgCode()) == 13:
                i_flav = 2 #numuCC
            elif abs(event.MCTrack[1].GetPdgCode()) == 15:
                is1Mu = False
                for j_track in range(2, len(event.MCTrack)):
                    if event.MCTrack[j_track].GetMotherId() == 1 and abs(event.MCTrack[j_track].GetPdgCode()) == 13:
                        is1Mu = True
                        break
                if is1Mu:
                    i_flav = 4 #nutauCC1mu
                else:
                    i_flav = 3 #nutauCC0mu            

        if BDT_cut:

            if ch.GetName()=='cbmsim':
                hcalTools.setsimulation(True) # Set simulation bool for hcalTools.muAna as well
                hcalTools.eventHasMuon=hcalTools.OutgoingMuon(event, 'neutrino')
                hcalTools.NeutrinoIntType=i_flav
            elif ch.GetName()=='rawConv':
                # hcalTools.eventHasMuon=np.nan
                runNr = event.EventHeader.GetRunId() 
                hcalTools.setsimulation(False, runNr) # For real data, ensure the correct timewalk and time alignment parameters are being used
            hcalTools.EventNumber=i_event 

            """
            Determines all the features needed for the BDT (list at: hcalTools.model.feature_names_in)
            To determine these features hcalTools is in the sndsw/python/ directory and utilises the barycentre and lambda methods in
            sndsw/python/AnalysisFunctions.py. These should eventually be incorporated into the MuFilter class. 
            """

            bdt_probs = hcalTools.zeromuBDTcut(event, muFilterDet,scifiDet)
            prob_noMuon, prob_Muon = bdt_probs[0,0], bdt_probs[0,1]
            
            BDT_weights[i_flav].Fill(prob_noMuon)
            BDT_weights[-1].Fill(prob_noMuon)
            
            if prob_noMuon <= hcalTools.BDTthreshold:  # If BDT assigns less than opt. thres. probability that event is signal (no final state muon) then ditch it
                print(f'BDT cut event {i_event}')
                continue

        weight = 1

        if isNeutralHad:
            weight = event.per_invfb_weight*LUMI

        elif isNuMC:
            weight = LUMI/(N_MC_FILES*100.)

        weight /= preselection

        if isNeutralHad:
            totWeight += weight
        
        n[0] += 1
        passCuts = True

        selHits = selectHits(event, (isNuMC or isNeutralHad))
        
        # Temporary fix! 
        if len(selHits)==1: continue

        slopev, slopeh, reducedchi2v, reducedchi2h, reducedchi2both = getSciFiAngle(selHits)
  
        h_min_chi2[i_flav].Fill(min(reducedchi2v, reducedchi2h), weight)
        h_log_min_chi2[i_flav].Fill(ROOT.TMath.Log(min(reducedchi2v, reducedchi2h)), weight)

        dens_sel, dens_sel2 = getSumDensity(selHits, return_2ndhighest = True)
        h_hit_density_sel_precut[i_flav].Fill(dens_sel, weight)
        h_hit_density2_sel_precut[i_flav].Fill(dens_sel, weight)

        if dens_sel2 < 20:
            continue

        h_hit_density_sel_after_dens2[i_flav].Fill(dens_sel, weight)

        if dens_sel < 2000:
            continue

        h_hit_density_sel_after_dens[i_flav].Fill(dens_sel, weight)
        
        if not (isNuMC or isNeutralHad):
            selected_list.append("SELECTED {} {} {} {}".format(i_event, event.EventHeader.GetRunId(), event.EventHeader.GetEventNumber(), dens_sel))

        h_SciFiAngle[i_flav].Fill(slopev, slopeh, weight)
        h_SciFiAngle_v_chi2[i_flav].Fill(slopev, reducedchi2v, weight)
        h_SciFiAngle_h_chi2[i_flav].Fill(slopeh, reducedchi2h, weight)

        theta = (ROOT.TMath.ATan(slopev)**2 + ROOT.TMath.ATan(slopeh)**2)**0.5

        h_theta[i_flav].Fill(theta, weight)
        
        h_n_hits_sel[i_flav].Fill(len(selHits), weight)
    
        n_hits = 0
        for hit in event.Digi_ScifiHits:
            if hit.isValid():
                n_hits += 1
        h_n_hits[i_flav].Fill(n_hits, weight)
    
        dens = getSumDensity(event.Digi_ScifiHits)[0]
        h_hit_density[i_flav].Fill(dens, weight)
        

        h_hit_density_sel[i_flav].Fill(dens_sel, weight)
        h_hit_density2_sel[i_flav].Fill(dens_sel2, weight)

        h_log_hit_density_sel[i_flav].Fill(ROOT.TMath.Log(dens_sel), weight)

        
        h_theta_density[i_flav].Fill(theta, dens_sel, weight)

        try:
            tdiff_list.append(tdiff(selHits, event.Digi_MuFilterHits, (isNuMC or isNeutralHad)))
        except ZeroDivisionError:
            tdiff_list.append(-999)
        max_slope_list.append(max(abs(slopev), abs(slopeh)))

    for i_sel, sel in enumerate(selected_list):
        print(sel, tdiff_list[i_sel], max_slope_list[i_sel])

    if isNeutralHad:
        for i_hist in (h_n_hits, h_n_hits_sel, h_hit_density, h_hit_density_sel, h_SciFiAngle, h_SciFiAngle_v_chi2, h_SciFiAngle_h_chi2, h_theta, h_theta_density, h_min_chi2, h_log_hit_density_sel, h_log_min_chi2, h_hit_density_sel_precut, h_hit_density2_sel, h_hit_density2_sel_precut, h_hit_density_sel_after_dens2, h_hit_density_sel_after_dens):
            for j_hist in i_hist:
                j_hist.Scale(N_HAD/totWeight)
        
    return (BDT_weights, h_n_hits, h_n_hits_sel, h_hit_density, h_hit_density_sel, h_SciFiAngle, h_SciFiAngle_v_chi2, h_SciFiAngle_h_chi2, h_theta, h_theta_density, h_min_chi2, h_log_hit_density_sel, h_log_min_chi2, h_hit_density_sel_precut, h_hit_density2_sel, h_hit_density2_sel_precut, h_hit_density_sel_after_dens2, h_hit_density_sel_after_dens)

out_file_name = f"checkDataCuts_BDTcut{options.BDTcut}.root"
out_file = ROOT.TFile(out_file_name, "recreate")

print(f'Making plots with data:\n')
plots_data  = makePlots(ch, options.BDTcut)
print(f'Making plots with neutrino MC:\n')
plots_MC = makePlots(chMC, options.BDTcut, isNuMC = True, name = "_MC")
print(f'Making plots with neutral hadron MC:\n')
plots_hadMC = makePlots(chNeutral, options.BDTcut, isNeutralHad = True, name = "_hadMC", preselection = 1.0)

out_file.Write()
out_file.Close()
print(f'File for BDT_cut {options.BDTcut} written to {out_file_name}')

def testing(ch, i):
    ch.GetEvent(i)
    hcalTools.EventNumber=i
    if ch.GetName()=='cbmsim':
        hcalTools.setsimulation(True)
        hcalTools.eventHasMuon=hcalTools.OutgoingMuon(ch)
    elif ch.GetName()=='rawConv':
        runNr = ch.EventHeader.GetRunId()
        hcalTools.setsimulation(False, runNr)

    data=hcalTools.BDT_cut(ch, muFilterDet, scifiDet)
    return data

c = []
def drawDataMC(data, MC):
    c.append(ROOT.TCanvas())
    c[-1].Divide(2, 1)
    c[-1].cd(1)
    data.Draw("COLZ")
    ROOT.gPad.Update()
    c[-1].cd(2)
    MC.Draw("COLZ")
    ROOT.gPad.Update()

#drawDataMC(h_SciFiAngle, h_SciFiAngle_MC)
#drawDataMC(h_SciFiAngle_v_chi2, h_SciFiAngle_v_chi2_MC)
#drawDataMC(h_SciFiAngle_h_chi2, h_SciFiAngle_v_chi2_MC)
#
#c.append(ROOT.TCanvas())
#h_n_hits.Draw("PE")
#if N_MC_FILES:
#    h_n_hits_MC.SetLineColor(ROOT.kRed)
#    h_n_hits_MC.Draw("SAME")
#    h_n_hits_sel_MC.SetLineColor(ROOT.kRed+2)
#    h_n_hits_sel_MC.Draw("SAME")
#ROOT.gPad.Update()
#
#c.append(ROOT.TCanvas())
#h_hit_density.Draw("PE")
#if N_MC_FILES:
#    h_hit_density_MC.SetLineColor(ROOT.kRed)
#    h_hit_density_MC.Draw("SAME")
#    h_hit_density_sel_MC.SetLineColor(ROOT.kRed+2)
#    h_hit_density_sel_MC.Draw("SAME")
#ROOT.gPad.Update()
#
#c.append(ROOT.TCanvas())
#h_theta.Draw("PE")
# if N_MC_FILES:
#    h_theta_MC.SetLineColor(ROOT.kRed)
#    h_theta_MC.Draw("SAME")
# ROOT.gPad.Update()
#
#c.append(ROOT.TCanvas())
#h_min_chi2.Draw("PE")
#if N_MC_FILES:
#    h_min_chi2_MC.SetLineColor(ROOT.kRed)
#    h_min_chi2_MC.Draw("SAME")
#ROOT.gPad.Update()
#
#for h in h_summary:
#    h.Write()

#input()
