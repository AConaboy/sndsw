import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import math
import json
from matplotlib.patches import Polygon
from skspatial.objects import Line, Points
from skspatial.plotting import plot_3d
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection
from mpl_toolkits.mplot3d import Axes3D
from scipy.stats import norm, skewnorm
from scipy.stats import poisson
from scipy.stats import chi2
from ipywidgets import interact, fixed
from collections import namedtuple
from collections import defaultdict
import pickle
import argparse, os, sys
import csv
import ROOT
pd.options.mode.copy_on_write = True

from snd_simulation_init import initialise_data
from snd_simulation_init import barycenter_calc
from snd_simulation_init import build_track
from snd_simulation_init import muon_reconstruction
from snd_simulation_init import FiducialCutProcessor
from snd_simulation_init import Tracking

parser = argparse.ArgumentParser(description="snd simulation")
parser.add_argument('--initdata', action="store_true", help='initiate data')
parser.add_argument('--track_reconstruction', action="store_true", help='track reconstruction')
parser.add_argument('--debug', action="store_true", help='debugging')
args = parser.parse_args()  
eospath = '/eos/user/t/tismith/SWAN_projects/genie_ana_output/'     

#get bar positions in US
with open(f'{eospath}BarPositions.json', 'r') as jf:
    barpositions = json.load(jf)
    barpositions = {int(k):v for k,v in barpositions.items()} # make sure the keys are integers not strings

#get z-positions of planes in DS
with open(f'{eospath}zPositions.data', 'rb') as f:
    zPos=pickle.load(f)
    
#c_scint for bars in US
#SiPMs 0-7 are on the left, 8-15 on the right
#dictionary: {detID_SiPM: [cscint, error]}
with open(f'{eospath}run005408_cscint_corrected.data', 'rb') as f:
    cscint=pickle.load(f)

fiducial_processor = FiducialCutProcessor(70, 105, 10, 50, 200, 1200, 300, 128*12-200, True, True)
setup_data = initialise_data(energy_loss_threshold=0.001, time_cut=31.25, attenuation_length=350)
bary_setup = barycenter_calc(energy_loss_threshold=0.001, QDCbar_threshold=0.002)
muon_recon = muon_reconstruction()
track = build_track(energy_loss_threshold=0.001)
tracking = Tracking()
tracking.Init()

background = True
amount_sel = 10
sel_cuts = True
#1: interaction in 5th target wall, 2: at least 2 consecutive scifi planes are hit, 3: if DS hits, then all US planes must be hit
#4: event has one reconstructed DS track -> in MC: clustering algorithm produces clusters and fit is possible
#5: latest DS hit time > earliest scifi hit time, 6: muon track intersects first scifi plane >5cm away from detector edge
#7: sum of min(DOCA) (=distance of closest approach) of track to scifi hits is <3cm in horizontal and vertical direction per station(plane in scifi)
#8: more than 35 scifi hits, 9: US total QDC larger than 700, 10: number of DS hits per projection >10
veto = -1

if args.initdata:
    signal = False
    muonDIS = True
    pm = False
    nh = False
    if signal == True:
        #initialise data numu
        #get HCAL data 
        dtype_options = {"Process": str}
        df = pd.read_csv(f"{eospath}data_us_10.csv", dtype=dtype_options).drop(columns = ["Unnamed: 0"]) #hcal data
        df = df.rename(columns={'FLUKA_weight': 'w'})
        a = df["Eventnumber"].nunique()
        print(f"Amount of events pre cuts: {a}")

        #get scifi data
        df_scifi = pd.read_csv(f"{eospath}data_scifi_10.csv", dtype=dtype_options).drop(columns = ["Unnamed: 0"]) #scifi data
        df_scifi["StationNR"] = df_scifi["detectorID"].apply(lambda x: math.floor(x / 1000000))
        a = df["Eventnumber"].nunique()
        cut_flow = []
        cut_flow.append(('Initial', a))

        #apply shower wall algo
        df_scifi["wall"] = None
        df_scifi = setup_data.shower_wall_algo(df_scifi)
        scifi_interaction_walls = df_scifi.groupby("Eventnumber")["wall"].first().reset_index()
        df = df.merge(scifi_interaction_walls, on="Eventnumber", how="left", suffixes=("", "_updated"))

        #apply data cuts
        #1cc = 6.25ns
        df_muon, df_shower, df, cut_flow = setup_data.init(df, df_scifi, fiducial_cuts=True, plot_fiducial=False, 
                                                        plot_time=False, process="numu", veto=veto, 
                                                        selection_cuts=sel_cuts, amount_selection_cuts=amount_sel, anti_cut=False, 
                                                        momentum_cut=False, slope_cut=False, cut_flow=cut_flow) 

        print(cut_flow)
        # filename = "10"
        # df.to_csv(f"{eospath}data_us_{filename}_fiducial_selCuts{amount_sel}.csv")
        # print(f'Numu csv written to: {eospath}data_us_{filename}_fiducial_selCuts{amount_sel}.csv')
        # filename = "10"
        # df.to_csv(f"{eospath}data_us_{filename}_fiducial.csv")
        # print(f'Numu csv written to: {eospath}data_us_{filename}_fiducial.csv')

    if background == True:
        if muonDIS == True:
            #get data for muon DIS
            #get HCAL data 
            df_muonDIS = pd.read_csv(f"{eospath}MuonDIS_all_weights.csv").drop(columns=["Unnamed: 0"])
            # df_muonDIS = pd.read_csv("MuonDIS_all.csv").drop(columns=["Unnamed: 0"])
            a = df_muonDIS["Eventnumber"].nunique()
            print(f"Amount of events pre cuts: {a}")
            cut_flow = []
            cut_flow.append(("Initial", a))

            #get scifi data
            df_muonDIS_scifi = pd.read_csv(f"{eospath}MuonDIS_all_scifi_weights.csv").drop(columns=["Unnamed: 0"])
            # df_muonDIS_scifi = pd.read_csv("MuonDIS_all_scifi.csv").drop(columns=["Unnamed: 0"])
            df_muonDIS_scifi["StationNR"] = df_muonDIS_scifi["detectorID"].apply(lambda x: math.floor(x / 1000000))

            #apply shower wall algo
            df_muonDIS_scifi["wall"] = None
            df_muonDIS_scifi = setup_data.shower_wall_algo(df_muonDIS_scifi)
            scifi_interaction_walls = df_muonDIS_scifi.groupby("Eventnumber")["wall"].first().reset_index()
            df_muonDIS = df_muonDIS.merge(scifi_interaction_walls, on="Eventnumber", how="left", suffixes=("", "_updated"))

            #apply data cuts
            df_muonDIS_muon, df_muonDIS_shower, df_muonDIS, cut_flow = setup_data.init(df_muonDIS, df_muonDIS_scifi, 
                                                                                    fiducial_cuts=True, plot_fiducial=False, 
                                                                                    plot_time=False, process="muonDIS", 
                                                                                    veto=veto, selection_cuts=sel_cuts,
                                                                                    amount_selection_cuts=amount_sel, anti_cut=False, 
                                                                                    momentum_cut=False, slope_cut=False, 
                                                                                    cut_flow=cut_flow)

            df_muonDIS.to_csv(f"{eospath}MuonDIS_fiducial_selCuts{amount_sel}.csv")
            print(f'MuonDIS csv written to: {eospath}MuonDIS_fiducial_selCuts{amount_sel}.csv')
            # df_muonDIS.to_csv(f"{eospath}MuonDIS_fiducial.csv")
            # print(f'MuonDIS csv written to: {eospath}MuonDIS_fiducial.csv')

        if nh == True:
            #get data for neutral hadrons
            #read HCAL data 
            df_nh = pd.read_csv(f"{eospath}NeutralHadrons_bigdata.csv").drop(columns=["Unnamed: 0"])
            a = df_nh["Eventnumber"].nunique()
            print(f"Amount of events pre cuts: {a}")
            cut_flow = []
            cut_flow.append(('Initial', a))

            #read scifi data
            df_nh_scifi = pd.read_csv(f"{eospath}NeutralHadrons_scifi_bigdata.csv").drop(columns=["Unnamed: 0"])
            df_nh_scifi["StationNR"] = df_nh_scifi["detectorID"].apply(lambda x: math.floor(x / 1000000))

            #apply shower wall algo
            df_nh_scifi["wall"] = None
            df_nh_scifi = setup_data.shower_wall_algo(df_nh_scifi)
            scifi_interaction_walls = df_nh_scifi.groupby("Eventnumber")["wall"].first().reset_index()
            df_nh = df_nh.merge(scifi_interaction_walls, on="Eventnumber", how="left", suffixes=("", "_updated"))

            #data cuts
            df_nh_muon, df_nh_shower, df_nh, cut_flow = setup_data.init(df_nh, df_nh_scifi, fiducial_cuts=True, 
                                                                        plot_fiducial=False, plot_time=False, process="nh", 
                                                                        veto=veto, selection_cuts=sel_cuts, 
                                                                        amount_selection_cuts=amount_sel,
                                                                        anti_cut=False, momentum_cut=False, slope_cut=False, 
                                                                        cut_flow=cut_flow)

            df_nh.to_csv(f"{eospath}nh_fiducial_bigdata_selCuts{amount_sel}.csv")
            print(f'NH csv written to: {eospath}nh_fiducial_bigdata_selCuts{amount_sel}.csv')
            # df_nh.to_csv(f"{eospath}nh_fiducial_cutsapplied_bigdata.csv")
            # print(f'NH csv written to: {eospath}nh_fiducial_cutsapplied_bigdata.csv')

        if pm == True:
            #get passing muon data
            #read HCAL data
            df_pm = pd.read_csv(f"{eospath}PassingMuons_large.csv").drop(columns=["Unnamed: 0"])
            #this csv file contains the requirement that a muon be in the veto
            # df_pm = pd.read_csv("PassingMuons_temp.csv").drop(columns=["Unnamed: 0"])
            a = df_pm["Eventnumber"].nunique()
            print(f"Amount of events pre cuts: {a}")
            cut_flow = []
            cut_flow.append(('Initial', a))

            #read scifi data
            df_pm_scifi = pd.read_csv(f"{eospath}PassingMuons_scifi_large.csv").drop(columns=["Unnamed: 0"])
            df_pm_scifi["StationNR"] = df_pm_scifi["detectorID"].apply(lambda x: math.floor(x / 1000000))

            #apply shower wall algo
            df_pm_scifi["wall"] = None
            df_pm_scifi = setup_data.shower_wall_algo(df_pm_scifi)
            scifi_interaction_walls = df_pm_scifi.groupby("Eventnumber")["wall"].first().reset_index()
            df_pm = df_pm.merge(scifi_interaction_walls, on="Eventnumber", how="left", suffixes=("", "_updated"))

            #apply data cuts
            df_pm_muon, df_pm_shower, df_pm, cut_flow = setup_data.init(df_pm, df_pm_scifi, fiducial_cuts=True, 
                                                                        plot_fiducial=False, plot_time=False, process="pm", 
                                                                        veto=veto, selection_cuts=sel_cuts, 
                                                                        amount_selection_cuts=amount_sel,
                                                                        anti_cut=False, momentum_cut=False, slope_cut=False, 
                                                                        cut_flow=cut_flow)

            df_pm.to_csv(f"{eospath}pm_fiducial_selCuts{amount_sel}.csv")
            print(f'pm csv written to: {eospath}pm_fiducial_selCuts{amount_sel}.csv')
            # df_pm.to_csv(f"{eospath}pm_fiducial.csv")
            # print(f'pm csv written to: {eospath}pm_fiducial.csv')

if args.track_reconstruction:
    dtype_options = {"Process": str}
    # df = pd.read_csv(f"{eospath}nh_fiducial.csv", dtype=dtype_options).drop(columns = ["Unnamed: 0"])
    # df = pd.read_csv(f"{eospath}MuonDIS_fiducial.csv", dtype=dtype_options).drop(columns = ["Unnamed: 0"])
    # df = pd.read_csv(f"{eospath}pm_fiducial.csv", dtype=dtype_options).drop(columns = ["Unnamed: 0"])
    # df = pd.read_csv(f"{eospath}data_us_10_cutsapplied.csv", dtype=dtype_options).drop(columns = ["Unnamed: 0"])
    df = pd.read_csv(f"{eospath}data_us_10_fiducial.csv", dtype=dtype_options).drop(columns = ["Unnamed: 0"])
    df_DS = df[df["StationNR"]==3]
    tracking = Tracking()
    tracking.Init()
    hcal_grouped = df_DS.groupby("Eventnumber")
    a = df_DS["Eventnumber"].nunique()
    i = 0
    reduced_chi = []
    event_list = []
    pdg_list = []
    for event, group in hcal_grouped:
        muon_points, muon_track, clusters, polar_angle = muon_recon.reconstruct_muon(group)
        trackcandidates = tracking.DStrack(event, clusters)
        
        if len(trackcandidates) > 0:
            trackcandidate = trackcandidates[0]
        elif len(trackcandidates) == 0:
            continue
        hitlist = []
        for key, value in trackcandidate.items():
            hitlist.append(value)

        if len(hitlist) == 0:
            print("No hits in hitlist, skipping this event.")
            continue
        elif len(hitlist) > 0:
            rc = tracking.fitTrack(event, hitlist)
            if rc == -1:
                i += 1
            
            # fitStatus = rc.getFitStatus()
            # print("Chi-squared / Ndf:", fitStatus.getChi2() / (fitStatus.getNdf() + 1E-15))
            # print("Fit Converged:", fitStatus.isFitConverged())
            # print("Number of Degrees of Freedom:", fitStatus.getNdf())
            # reduced_chi.append(fitStatus.getChi2() / (fitStatus.getNdf() + 1E-15))
            # event_list.append(event)

            fittedState = rc.getFittedState()
            # covarianceMatrix = fittedState.get6DCov()
            # covarianceMatrix.Print()
            z = 379.70117959012583 # first US plane
            point = ROOT.TVector3(-40, 40, z)
            pos = fittedState.getPos()
            direction = fittedState.getDir()
            state = fittedState.get6DState()
            # extrap = fittedState.extrapolateToPoint(point)
            # pdg = fittedState.getPDG()
            print(event)
            # extrap.getPos().Print()
            pos.Print()
            direction.Print()
            state.Print()
            # print(pdg)
            # print(type(extrap))
            

    # print(f"Total amount of events: {a}")
    # print(f"Events with converged tracks: {a-i}")

    # track_data = pd.DataFrame()
    # track_data["Eventnumber"] = event_list
    # track_data["reduced_chi"] = reduced_chi

    # track_data.to_csv(f"{eospath}track_info_nh.csv")
    # print(f'Track info csv written to: {eospath}track_info_nh.csv')

if args.debug:

    #get HCAL data 
    dtype_options = {"Process": str}
    # df = pd.read_csv(f"{eospath}NeutralHadrons_bigdata.csv", dtype=dtype_options).drop(columns = ["Unnamed: 0"]) #hcal data
    # df_scifi = pd.read_csv(f"{eospath}NeutralHadrons_scifi_bigdata.csv", dtype=dtype_options).drop(columns = ["Unnamed: 0"])
    # df_scifi["StationNR"] = df_scifi["detectorID"].apply(lambda x: math.floor(x / 1000000))

    # cut_flow = []
    # a = df["Eventnumber"].nunique()
    # print(f"Amount of events pre cuts: {a}")
    # cut_flow.append(("Initial", a))

    # #apply shower wall algo
    # df_scifi["wall"] = None
    # df_scifi = setup_data.shower_wall_algo(df_scifi)
    # scifi_interaction_walls = df_scifi.groupby("Eventnumber")["wall"].first().reset_index()
    # df = df.merge(scifi_interaction_walls, on="Eventnumber", how="left", suffixes=("", "_updated"))

    # veto = -1
    # amount_sel = 10
    # sel_cuts = True
    # fiducial = True
    # #apply data cuts
    # #1cc = 6.25ns
    # df_muon, df_shower, df, cut_flow = setup_data.init(df, df_scifi, fiducial_cuts=fiducial, plot_fiducial=False, 
    #                                                 plot_time=False, process="pm", veto=veto, 
    #                                                 selection_cuts=sel_cuts, amount_selection_cuts=amount_sel, anti_cut=False, 
    #                                                 momentum_cut=False, slope_cut=False, cut_flow=cut_flow) 

    #add track info 
    # df["track_pos"] = None
    # df["track_dir"] = None
    # df["track_chi2"] = None
    # df["track_dof"] = None
    # df["track_covM"] = None
    
    # filename = "data_us_10_fiducial.csv"
    filename = "nh_fiducial_bigdata.csv"
    df = pd.read_csv(f"{eospath}{filename}", dtype=dtype_options).drop(columns = ["Unnamed: 0"]) #hcal data
    i = 0
    df_DS = df[df["StationNR"]==3]
    hcal_grouped = df_DS.groupby("Eventnumber")
    for event, group in hcal_grouped:
        muon_points, muon_track, clusters, polar_angle = muon_recon.reconstruct_muon(group)
        trackcandidates = tracking.DStrack(event, clusters)
        
        if len(trackcandidates) > 0:
            trackcandidate = trackcandidates[0]
        elif len(trackcandidates) == 0:
            continue
        hitlist = []
        for key, value in trackcandidate.items():
            hitlist.append(value)

        if len(hitlist) == 0:
            continue
        elif len(hitlist) > 0:
            rc = tracking.fitTrack(event, hitlist)
            if rc == -1:
                i += 1
                continue

            fittedState = rc.getFittedState()
            fitStatus = rc.getFitStatus()
            chi2 = fitStatus.getChi2()
            dof = fitStatus.getNdf() + 1E-15
            reduced_chi2 = chi2 / dof
            covarianceMatrix = fittedState.get6DCov()
            pos = fittedState.getPos()
            direction = fittedState.getDir()

            pos = [pos.X(), pos.Y(), pos.Z()]
            direction = [direction.X(), direction.Y(), direction.Z()]
            covarianceMatrix_flat = [covarianceMatrix(i, j) for i in range(6) for j in range(6)]

            # df.loc[df["Eventnumber"] == event, "track_pos"] = pos 
            df.loc[df["Eventnumber"] == event, ["track_posX", "track_posY", "track_posZ"]] = [pos[0], pos[1], pos[2]]
            df.loc[df["Eventnumber"] == event, ["track_dirX", "track_dirY", "track_dirZ"]] = [direction[0], direction[1], direction[2]]
            df.loc[df["Eventnumber"] == event, "track_chi2"] = chi2
            df.loc[df["Eventnumber"] == event, "track_dof"] = dof
            # df.loc[df["Eventnumber"] == event, "track_covM"] = covarianceMatrix_flat

    df.to_csv(f"{eospath}track_info_nh.csv")
    # print(cut_flow)
    # filename = "10"
    # df.to_csv(f"{eospath}nh_fiducial_selCuts{amount_sel}_NEW.csv")
    # print(f'csv written to: {eospath}nh_fiducial_selCuts{amount_sel}_NEW.csv')
    # with open(f'{eospath}nh_cutflow_fiducial_selCuts{amount_sel}_NEW.csv', 'w', newline='') as f:
    #     writer = csv.writer(f)
    #     writer.writerows(cut_flow)
