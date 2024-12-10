import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import math
import json
from array import array
from matplotlib.patches import Polygon
from skspatial.objects import Line, Points
from skspatial.plotting import plot_3d
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection
from mpl_toolkits.mplot3d import Axes3D
from scipy.stats import norm, skewnorm
from scipy.stats import poisson
from scipy.stats import linregress
from scipy.stats import chi2
from ipywidgets import interact, fixed
from collections import namedtuple
from collections import defaultdict
import pickle
import argparse, os, sys
import ROOT 
import SndlhcGeo
pd.options.mode.copy_on_write = True

class initialise_data:
    def __init__(self, energy_loss_threshold, time_cut, attenuation_length):
        self.energy_loss_threshold = energy_loss_threshold
        self.time_cut = time_cut
        self.attenuation_length = attenuation_length

    def init(self, df, df_scifi, fiducial_cuts, plot_fiducial, plot_time, process, veto, 
             selection_cuts, amount_selection_cuts, anti_cut, momentum_cut, slope_cut, cut_flow):
        fiducial_processor = FiducialCutProcessor(70, 105, 10, 50, 200, 1200, 300, 128*12-200, True, True)
        muon_recon = muon_reconstruction()
        
        #get detector dimensions
        hcal_xmin, hcal_xmax = df["coordX"].min(), df["coordX"].max()
        hcal_ymin, hcal_ymax = df["coordY"].min(), df["coordY"].max()
        scifi_xmin, scifi_xmax = df_scifi["coordX"].min(), df_scifi["coordX"].max()
        scifi_ymin, scifi_ymax = df_scifi["coordY"].min(), df_scifi["coordY"].max()
        dim_detector = [hcal_xmin, hcal_xmax, hcal_ymin, hcal_ymax, scifi_xmin, scifi_xmax, scifi_ymin, scifi_ymax]
        
        if process == "numu":
            # Drop outliers in time
            df = self.cut_outliers_t(df, plot=plot_time)
            cut_flow.append(('time interval ' f"[0, ({self.time_cut}ns)]", df["Eventnumber"].nunique()))
    #             cut_flow.append(("Sim cut 1", df["Eventnumber"].nunique()))

        # Drop heavy particles
        df = df[df["pdg_code"] < 1e5]
        cut_flow.append(('heavy elements', df["Eventnumber"].nunique()))
#             cut_flow.append(('Sim cut 2', df["Eventnumber"].nunique()))
            
        if process == "numu":
            # Remove all events without muon/anti-muon
            tot = df["Eventnumber"].nunique()
            muon_mask = df.groupby("Eventnumber")["pdg_code"].transform(lambda x: 13 in x.values or -13 in x.values)
            df = df[muon_mask]
            b = df["Eventnumber"].nunique()
            a = tot - b
            cut_flow.append(('muon presence', b))
#             cut_flow.append(('Sim cut 3', b))
            print("event contains muon: " + str(b) + ", Ratio NC/CC = " + str(round(a/b, 4) * 100) + "%")

        # Smear energy using the attenuation length
        df["Energy_loss"] *= np.exp(-df["coordX"] / self.attenuation_length)

        # Drop MC points that don't pass the energy threshold
        df = df[df["Energy_loss"] > self.energy_loss_threshold]
        cut_flow.append(('MC point energy threshold (1MeV)', df["Eventnumber"].nunique()))
#         cut_flow.append(('Sim cut 4', df["Eventnumber"].nunique()))

        # Translate detectorID
        df["StationNR"] = df["detectorID"] // 10000
        df["PlaneNR"] = (df["detectorID"] % 10000) // 1000
        df["Bar"] = df["detectorID"] % 1000

        # Apply fiducial cuts
        if fiducial_cuts:
            df_preFiducial = df
            df, df_scifi, pass_hcal, pass_scifi = fiducial_processor.apply_cuts(df, df_scifi, plot_fiducial)
            if pass_scifi > pass_hcal:
                cut_flow.append(('fiducial cuts: SciFi', pass_scifi))
                cut_flow.append(('fiducial cuts: DS', pass_hcal))
#                 cut_flow.append(('Selection cut A', pass_scifi))
#                 cut_flow.append(('Selection cut B', pass_hcal))
            else:
                cut_flow.append(('fiducial cuts: DS', pass_hcal))
                cut_flow.append(('fiducial cuts: SciFi', pass_scifi))
#                 cut_flow.append(('Selection cut B', pass_hcal))
#                 cut_flow.append(('Selection cut A', pass_scifi))
        
        # Separate muon and shower
        df_muon = df[((df["pdg_code"] == 13) | (df["pdg_code"] == -13))]
        df_shower = df[(df["pdg_code"] != 13) & (df["pdg_code"] != -13)]
            
        #remove events if veto fires
        events_firstplane = df_scifi[df_scifi["StationNR"] == 1]["Eventnumber"].unique()
        events_secondplane = df_scifi[df_scifi["StationNR"] == 2]["Eventnumber"].unique()
        if veto >= 0:
            events_with_stationnr_1 = df[df["StationNR"] == 1]["Eventnumber"].unique()
            df['random_number'] = np.random.uniform(0, 1, len(df)) #random number for veto inefficiency
            condition = (~df['Eventnumber'].isin(events_with_stationnr_1)) & (df['random_number'] > 5e-5)
            df = df[condition]
            cut_flow.append(('veto fires', df["Eventnumber"].nunique()))
#             cut_flow.append(('Selection cut C', df["Eventnumber"].nunique()))
            if veto == 1 or veto == 2:
                df['random_number'] = np.random.uniform(0, 1, len(df)) #random number for scifi inefficiency
#                 df = df[~df["Eventnumber"].isin(events_firstplane)]
                condition = (~df['Eventnumber'].isin(events_firstplane)) & (df['random_number'] > 1.1e-4)
                df = df[condition]
                cut_flow.append(("scifi plane 1 fires", df["Eventnumber"].nunique()))
#                 cut_flow.append(("Selection cut D", df["Eventnumber"].nunique()))
            if veto == 2:
#                 df = df[~df["Eventnumber"].isin(events_secondplane)]
                df['random_number'] = np.random.uniform(0, 1, len(df)) #random number for scifi inefficiency
                condition = (~df['Eventnumber'].isin(events_secondplane)) & (df['random_number'] > 1.1e-4)
                df = df[condition]
                cut_flow.append(("scifi plane 2 fires", df["Eventnumber"].nunique()))
#                 cut_flow.append(("Selection cut E", df["Eventnumber"].nunique()))
        
        if selection_cuts == True:
            if amount_selection_cuts >= 1:
                #remove events with interaction vertex in 5th wall
                df = df[df["wall"] != 4]
                cut_flow.append(("interaction in 5th wall", df["Eventnumber"].nunique()))
        #             cut_flow.append(("Selection cut F", df["Eventnumber"].nunique()))

            #perform selection cuts
            df, df_scifi, cut_flow = self.selection_cuts(df, df_scifi, dim_detector, cut_flow, 
                                                         process=process, amount_selection_cuts=amount_selection_cuts)
    #         cut_flow.append(("selection cuts", df["Eventnumber"].nunique()))
    
        if anti_cut == True:
            hcal_grouped = df.groupby("Eventnumber")
            events_to_remove = []
            for event, group in hcal_grouped:
                muon_points, muon_track, clusters, polar_angle = muon_recon.reconstruct_muon(group)
                if polar_angle != None:
                    angle = math.degrees(polar_angle)
                    if angle <= 30:
                        events_to_remove.append(event)
            df = df[~df["Eventnumber"].isin(events_to_remove)]        
            cut_flow.append((r"anti cut: $\theta$ > 30", df["Eventnumber"].nunique()))

        # Separate muon and shower
        df_muon = df[((df["pdg_code"] == 13) | (df["pdg_code"] == -13))]
        df_shower = df[(df["pdg_code"] != 13) & (df["pdg_code"] != -13)]
        
        if momentum_cut or slope_cut:
            #perform cuts on muons
            df_muon, momentum_cut, slope_cut, cut_flow = self.muon_cuts(df_muon, cut_flow, 
                                                                        momentum_cut=momentum_cut, slope_cut=slope_cut)
            df_muon = df_muon.reset_index(drop=True)
            df = pd.concat([df_shower, df_muon])

        return df_muon, df_shower, df, cut_flow
    
    def plot_cutflow(self, cut_flow, title):
        cut_names, cut_values = zip(*cut_flow)
        plt.figure(figsize=(15, 10))
#         plt.bar(range(len(cut_names)), cut_values, tick_label=cut_names)
        bars = plt.barh(cut_names, cut_values)
        plt.gca().invert_yaxis()
        plt.xlabel('Number of Unique Events', fontsize=16)
        plt.title(f"{title}" " cut flow, events left: " f"{cut_values[len(cut_values)-2]}", fontsize=16)
        plt.xticks(rotation=45, ha='right', fontsize=16)
        plt.yticks(fontsize=16)
        plt.xlim(0, cut_values[0]*1.25)
        
#         x_min, x_max = plt.xlim()
#         x_center = (x_min + x_max) / 2
        for bar, value in zip(bars, cut_values):
                plt.text(cut_values[0]*1.15, bar.get_y() + bar.get_height() / 2,  # Center text horizontally and vertically
                         f"{value} events", fontsize=16, color="black", ha='center', va='center')
        plt.show()

    #cut outliers in time
    def cut_outliers_t(self, df, plot):

        df_valid = df[df["time"]<self.time_cut]
        
        print("time intervall before cut: [0, " + str(round(df["time"].max(), 3)) + "]\ntime intervall after cut: [0, " 
              + str(round(df_valid["time"].max(), 3)) + "]")
        
        if plot == True:
            #shifted time for shower and event, starting time of muon at t=1
            df_valid["time_shifted"] = None
            for e in df_valid["Eventnumber"].unique():
                temp = df_valid.loc[df_valid["Eventnumber"]==e]
                min_time = temp["time"].min()
                df_valid.loc[df_valid["Eventnumber"] == e, "time_shifted"] = temp["time"] - min_time + 1
    
            time_precut, time_postcut = df["time"], df_valid["time"]
            T, X = df_valid["time"], df_valid["coordX"]
            T_shifted = df_valid["time_shifted"]
            
            #plot x-t-histogram
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), sharey=True)

            h1 = ax1.hist2d(T_shifted, X, bins=300, norm=mpl.colors.LogNorm(), cmap=plt.cm.jet)
            cbar1 = fig.colorbar(h1[3], ax=ax1)
            cbar1.set_label('Counts', fontsize=14)
            ax1.set_xlabel('time (ns)', fontsize=14)
            ax1.set_ylabel('x-coordinate', fontsize=14)
            ax1.set_title("time distribution")
            ax1.set_xlim(0, self.time_cut)
            
#             time_lines = [1, 1.45, 1.9, 2.35, 2.8]
#             for time in time_lines:
#                 ax2.axvline(time, color="black", linewidth=1.25)
            h2 = ax2.hist2d(T_shifted, X, bins=500, norm=mpl.colors.LogNorm(), cmap=plt.cm.jet)
            cbar2 = fig.colorbar(h2[3], ax=ax2)
            cbar2.set_label('Counts', fontsize=14)
            ax2.set_xlim(0.5, 10)
            ax2.set_ylim(-55,0)
            ax2.set_xlabel('time (ns)', fontsize=14)
            ax2.set_title("time distribution with xlim")

            plt.tight_layout()
            plt.show()

            #time distributions
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4), sharey=True)

            #pre cut
            bins_precut = np.arange(0, 1e3, 0.5)
            ax1.hist(time_precut, bins=bins_precut)
            ax1.set_xlim(0, 1e3)
            ax1.set_ylabel("Counts", fontsize=14)
            ax1.set_xlabel("time (ns)", fontsize=14)
            ax1.set_title("time distribution pre cut", fontsize=14)

            #post cut
            bins_postcut = np.arange(0, self.time_cut, 0.5)
            ax2.hist(time_postcut, bins=bins_postcut)
            ax2.set_xlabel("time (ns)", fontsize=14)
            ax2.set_title("time distribution post cut", fontsize=14)

            # Display the plots
            plt.tight_layout()
            plt.show()

        return df_valid
    
    def read_quark_data(self, df, filename):
        # Read data from CSV file and drop the 'Unnamed: 0' column
        df_read_quark = pd.read_csv(str(filename)).drop(columns=["Unnamed: 0"], errors="ignore")

        # Transpose the DataFrame and reset index
        df_quark_transposed = df_read_quark.T.reset_index()

        # Rename columns to 'eventnumber' and 'particle_info_0', 'particle_info_1', ...
        df_quark_transposed.columns = ['Eventnumber'] + [f'particle_info_{i}' for i in 
                                                         range(len(df_quark_transposed.columns) - 1)]
        df_quark_transposed['Eventnumber'] = df_quark_transposed['Eventnumber'].astype(int)

        # Convert 'particle_info' columns from string to dictionary
        for col in df_quark_transposed.columns[1:]:
            df_quark_transposed[col] = df_quark_transposed[col].apply(eval)

        # Flatten the DataFrame
        flattened_data = [
            {**row[f'particle_info_{i}'], 'Eventnumber': row['Eventnumber']}
            for _, row in df_quark_transposed.iterrows()
            for i in range(len(df_quark_transposed.columns) - 1)]
        
        df_quark_info = pd.DataFrame(flattened_data)
        df_numu = df_quark_info
#         df_numu = df_quark_info.loc[(df_quark_info["pdg"]==14) | (df_quark_info["pdg"]==-14)]
        df_numu = df_numu[df_numu["Eventnumber"].isin(df["Eventnumber"].unique())]
        df_numu = df_numu.reset_index(drop=True)

        return df_numu
    
    def shower_wall_algo(self, df):
        k = 0.1 #parameter used as threshold
        df["wall"] = None
        
        for event in df["Eventnumber"].unique():
            df_event = df[df["Eventnumber"]==event]
            i = 1
            k_temp = 0
            energy_loss_plane = df_event.groupby("StationNR")["Energy_loss"].sum()
#             print(event, energy_loss_plane)
            while k_temp <= k and i <= 5:
                if i not in energy_loss_plane or i + 1 not in energy_loss_plane:
#                     if event == 2072:
#                         print(i, "hi")
                    i += 1
                else:
                    E_wall = energy_loss_plane[i]
                    E_wallnext = energy_loss_plane[i+1]
                    k_temp = (abs(E_wall-E_wallnext))/(E_wall)
#                     if event == 2072:
#                         print(E_wall, E_wallnext, k_temp)
                    if k_temp > k:
                        df.loc[(df["Eventnumber"]==event), "wall"] = i  
#                         if event == 2072:
#                             print(i, "wall saved")
                    i += 1
        return df
    
    #cuts on muon momentum and slope
    def muon_cuts(self, df_muon, cut_flow, momentum_cut, slope_cut):
        plot_dict = {"hist_values":[], "hist_label":[], "x_label":[], "y_label":[], "title":[]}
        event_amount = df_muon["Eventnumber"].nunique()
        
        if momentum_cut == True:
            #cuts on muon momenta
            muon_pz = df_muon["pz"]
            muon_pz_decays = df_muon.loc[(df_muon["TrackID"] > 1)&(df_muon["pz"] < 600), "pz"]
            momentum_cut = 25
            print(f"max muon momentum: {momentum_cut} MeV")

            #perform muon cuts based on momenta in z
            df_muon = df_muon.loc[(df_muon["pz"] > momentum_cut)&(df_muon["pz"]<600)]
            a = df_muon["Eventnumber"].nunique()
            cut_flow.append(("muon momentum cut: "+"$p_z = $"+f"{round(momentum_cut, 2)}MeV", a))
    #         cut_flow.append(("Sim cut 5", a))

        #cuts on muon slope 
        slope = np.degrees(np.arctan2(df_muon["px"], df_muon["pz"] ))
        cut_slope_degrees = 10

        if slope_cut == True:
            mask = np.abs(slope) <= cut_slope_degrees
            df_muon = df_muon[mask] 

            a = df_muon["Eventnumber"].nunique()
            cut_flow.append((f"muon slope cut: {cut_slope_degrees}°", a))
    #         cut_flow.append(("Sim cut 6", a))
        
        print("Amount of events pre muon cuts: " + str(event_amount))
        print("Amount of events post muon cuts: " + str(df_muon["Eventnumber"].nunique()))

        return df_muon, momentum_cut, cut_slope_degrees, cut_flow
    
    def plot_muon_cuts(self, hist_values, hist_label, x_label, y_label, title, i, momentum_cut):
        fig = plt.figure(figsize=(8, 6))
        ax = fig.add_subplot(111)
        plt.hist(hist_values, bins = 100, label=hist_label)
        if i == 0:
            plt.axvline(momentum_cut, color="black", linestyle="--", 
                        label=r"proposed cut at $p_Z=$" + str(round(momentum_cut, 2)) + "MeV")
        plt.title(title, fontsize=14)
        plt.ylabel(y_label, fontsize=14)
        plt.xlabel(x_label, fontsize=14)
        plt.legend(fontsize=14)
        plt.show()
    
    #selection cuts as found in publication
    def selection_cuts(self, df, df_scifi, dim_detector, cut_flow, process, amount_selection_cuts):
        muon_recon = muon_reconstruction()
        
        hcal_grouped = df.groupby("Eventnumber")
#         df_scifi.rename(columns={"PlaneNR": "StationNR"}, inplace=True)
        scifi_grouped = df_scifi.groupby("Eventnumber")
        events_to_remove = []
        
        if amount_selection_cuts >= 2:
            #at least 2 consecutive scifi planes are hit
            for (event), group in scifi_grouped:
                i = 0
                for plane in group["StationNR"].unique():
                    group_plane = group[group["StationNR"]==plane]
                    if group_plane.empty:
                        i = 0
                    else:
                        i += 1 
                if i+1 < 2:
                    events_to_remove.append(event)

            temp = df[~df["Eventnumber"].isin(events_to_remove)]        
            cut_flow.append(("2 consecutive scifi planes hit", temp["Eventnumber"].nunique()))
    #         cut_flow.append(("Selection cut G", temp["Eventnumber"].nunique()))
        
        if amount_selection_cuts >= 3:
            #if there are DS hits, then all US planes must be hit
            for event, group in hcal_grouped:
                group_US = group[group["StationNR"] == 2]
                group_DS = group[group["StationNR"] == 3]
                if not group_DS.empty:
                    USplanes = group_US["PlaneNR"].nunique()
                    if not USplanes == 5 and event not in events_to_remove:
                        events_to_remove.append(event)
                elif event not in events_to_remove:
                    events_to_remove.append(event)

            temp = df[~df["Eventnumber"].isin(events_to_remove)]   
            cut_flow.append(("if DS hit: all US planes hit", temp["Eventnumber"].nunique()))
    #             cut_flow.append(("Selection cut H", temp["Eventnumber"].nunique()))
        
        if amount_selection_cuts >= 4:
            #event has one reconstructed DS track -> in MC: clustering algorithm produces clusters and fit is possible
            track = build_track(energy_loss_threshold=0.001)
            for event, group in hcal_grouped:
                muon_points, muon_track, clusters, angle = muon_recon.reconstruct_muon(group)
                if muon_track == None:
                    if event not in events_to_remove:
                        events_to_remove.append(event)

            temp = df[~df["Eventnumber"].isin(events_to_remove)]   
            cut_flow.append(("reconstructed muon track", temp["Eventnumber"].nunique()))
    #             cut_flow.append(("Selection cut J", temp["Eventnumber"].nunique()))
        
        if amount_selection_cuts >= 5:
            #latest DS hit time > earliest scifi hit time
            if process != "muonDIS":
                for event, group in scifi_grouped:
                    hcal_DS = df[(df["Eventnumber"]==event)&(df["StationNR"]==3)]
                    timeMIN_scifi = group["time"].min()
                    timeMAX_DS = hcal_DS["time"].max()
                    if timeMAX_DS <= timeMIN_scifi and event not in events_to_remove:
                        events_to_remove.append(event)

                temp = df[~df["Eventnumber"].isin(events_to_remove)]        
                cut_flow.append(("latest DS t > earliest scifi t", temp["Eventnumber"].nunique()))
        #             cut_flow.append(("Selection cut K", temp["Eventnumber"].nunique()))
        
        if amount_selection_cuts >= 6:
            #muon track intersects first scifi plane >5cm away from detector edge
            for event, group in hcal_grouped:
                muon_points, muon_track, clusters, angle = muon_recon.reconstruct_muon(group)
                if muon_track == None:
                    if event not in events_to_remove:
                        events_to_remove.append(event)
                else:
    #                     x_0, y_0, z_0 = muon_track.point
    #                     v_x, v_y, v_z = muon_track.direction
                    scifi_plane = zPos["Scifi"][10]
                    x_pos = muon_track['line_fit_x'](scifi_plane)
                    y_pos = muon_track['line_fit_y'](scifi_plane)
    #                     t = (scifi_plane - z_0) / v_z
    #                     x_pos = x_0 + t * v_x
    #                     y_pos = y_0 + t * v_y
                    x_diff, y_diff = True, True
                    x_diff1, x_diff2 = abs(x_pos-dim_detector[4]), abs(x_pos-dim_detector[5])
                    y_diff1, y_diff2 = abs(y_pos-dim_detector[6]), abs(y_pos-dim_detector[7])
                    if (x_diff1 < 5 or x_diff2 < 5 or y_diff1 < 5 or y_diff2 < 5) and event not in events_to_remove:
                        events_to_remove.append(event)

            temp = df[~df["Eventnumber"].isin(events_to_remove)]        
            cut_flow.append(("track intersects first scifi plane \n>5cm away from detector edge", temp["Eventnumber"].nunique()))
    #             cut_flow.append(("Selection cut L", temp["Eventnumber"].nunique()))
        
        if amount_selection_cuts >= 7:
            #sum of min(DOCA) (=distance of closest approach) of track to scifi hits is <3cm in horizontal and vertical
            #direction per station(plane in scifi)            
            for event, group in scifi_grouped:
                muon_points, muon_track, clusters, angle = muon_recon.reconstruct_muon(group)

                if muon_track is None:
                    if event not in events_to_remove:
                        events_to_remove.append(event)
                else:
    #                     x_0, y_0, z_0 = muon_track.point
    #                     v_x, v_y, v_z = muon_track.direction
                    station = 10
                    sum_per_station = []

                    while station < 52:
                        scifi_plane_z = zPos["Scifi"][station]
                        x_pos = muon_track['line_fit_x'](scifi_plane_z)
                        y_pos = muon_track['line_fit_y'](scifi_plane_z)
    #                         t = (scifi_plane_z - z_0) / v_z
    #                         x_pos = x_0 + t * v_x
    #                         y_pos = y_0 + t * v_y

                        # Filter scifi hits near the expected position (within some tolerance)
                        tolerance = 10
                        nearby_hits = df_scifi[
                            (np.abs(df_scifi['coordX'] - x_pos) < tolerance) & 
                            (np.abs(df_scifi['coordY'] - y_pos) < tolerance)
                        ]

                        if not nearby_hits.empty:
                            # Compute DOCA (distance of closest approach) in a vectorized way
                            distances = np.sqrt(
                                (nearby_hits['coordX'] - x_pos) ** 2 + 
                                (nearby_hits['coordY'] - y_pos) ** 2
                            )
                            min_distance = distances.min()
                            sum_per_station.append(min_distance)

                        station = station + 1 if station % 2 == 0 else station + 9

                    any_greater_than_3 = any(dist > 3 for dist in sum_per_station)
                    if any_greater_than_3 and event not in events_to_remove:
                        events_to_remove.append(event)

            temp = df[~df["Eventnumber"].isin(events_to_remove)]        
            cut_flow.append(("sum of min(DOCA) to scifi hits <3cm", temp["Eventnumber"].nunique()))
    #             cut_flow.append(("Selection cut M", temp["Eventnumber"].nunique()))
        
        if amount_selection_cuts >= 8:
            #more than 35 scifi hits
            for event, group in scifi_grouped:
                hits = len(group)
                if hits <= 35 and event not in events_to_remove:
                    events_to_remove.append(event)

            temp = df[~df["Eventnumber"].isin(events_to_remove)]        
            cut_flow.append((">35 scifi hits", temp["Eventnumber"].nunique()))
    #             cut_flow.append(("Selection cut N", temp["Eventnumber"].nunique()))

        if amount_selection_cuts >= 9:
            #US total QDC larger than 700
            #QDC calc: QDC calibration constant/number of SIPMs
            #QDC calibration constant: 25000
            for event, group in hcal_grouped:
                hcal_US = group[group["StationNR"]==2]
                hcal_US_grouped = hcal_US.groupby(["PlaneNR", "Bar"])
                sum_QDC_hcal = 0
                for (plane, bar), group2 in hcal_US_grouped:
                    sum_eloss_bar = group2["Energy_loss"].sum()
                    sum_QDC_bar = sum_eloss_bar*25*1000
                    sum_QDC_hcal += sum_QDC_bar
    #             print(sum_QDC_hcal)
                if sum_QDC_hcal < 700 and event not in events_to_remove:
                    events_to_remove.append(event)

            temp = df[~df["Eventnumber"].isin(events_to_remove)]        
            cut_flow.append(("sum US QDC >700", temp["Eventnumber"].nunique()))
    #             cut_flow.append(("Selection cut O", temp["Eventnumber"].nunique()))

        if amount_selection_cuts >= 10:
            #number of DS hits per projection >10 (I think just means event, see selection cuts scripts from Cristovao)
            for event, group in hcal_grouped:
                hcal_DS = group[group["StationNR"]==3]
                if len(hcal_DS) > 10 and event not in events_to_remove:
                    events_to_remove.append(event)

            temp = df[~df["Eventnumber"].isin(events_to_remove)]        
            cut_flow.append(("NR. of DS hits per event <10", temp["Eventnumber"].nunique()))
    #             cut_flow.append(("Selection cut P", temp["Eventnumber"].nunique()))

        df = df[~df["Eventnumber"].isin(events_to_remove)]
        df_scifi = df_scifi[~df_scifi["Eventnumber"].isin(events_to_remove)]
        
        return df, df_scifi, cut_flow            

#get bar positions in US
with open(f'/eos/user/t/tismith/SWAN_projects/genie_ana_output/BarPositions.json', 'r') as jf:
    barpositions = json.load(jf)
    barpositions = {int(k):v for k,v in barpositions.items()} # make sure the keys are integers not strings

#get z-positions of planes in DS
with open(f'/eos/user/t/tismith/SWAN_projects/genie_ana_output/zPositions.data', 'rb') as f:
    zPos=pickle.load(f)
    
#c_scint for bars in US
#SiPMs 0-7 are on the left, 8-15 on the right
#dictionary: {detID_SiPM: [cscint, error]}
with open(f'/eos/user/t/tismith/SWAN_projects/genie_ana_output/run005408_cscint_corrected.data', 'rb') as f:
    cscint=pickle.load(f)

#barycenter, residual and XY-plots in US and DS
class barycenter_calc:
    def __init__(self, energy_loss_threshold, QDCbar_threshold):
        self.threshold = energy_loss_threshold
        self.QDCbar_threshold = QDCbar_threshold

    def barycenter(self, df, p, process):
        c_scint = 14.5 #default constant scintillator velocity
        amount = 10 #set amount of bars in plane
        
        Result = namedtuple('Result', [
        'x_bary', 'y_bary', 'y_BarPositions', 
        'x_BarPositions', 'x_res', 'y_res', 'QDC_perBar', 'lambdaX', 'lambdaY',
        'mean_px', 'mean_pz', "interaction_wall", "weight", "lambdaY2"
        ]) #create tuple to return 
        
        df_bary = df.loc[(df["StationNR"] == 2) & (df["PlaneNR"] == p)]
        QDC_tot = df_bary["Energy_loss"].sum()  # entire energy loss in the plane
        y_min, y_max = df_bary["coordY"].mean(), df_bary["coordY"].mean()  # used to calculate lambda
        YMIN_global, YMAX_global = df_bary["coordY"].min(), df_bary["coordY"].max()
        
        #get interaction wall
        if process != "nh":
            unique_wall_values = df_bary["wall"].unique()
            unique_wall_values_list = unique_wall_values.tolist()
            if unique_wall_values_list: 
                unique_wall_value = unique_wall_values_list[0]
                interaction_wall = -1 if pd.isna(unique_wall_value) else int(unique_wall_value)
            else:
                interaction_wall = -1
        else:
            interaction_wall = -1
        
        results = {
        'x': [], #x barycenter
        'y': [], #y barycenter
        'y_BarPositions': [], #middle position of bars in y
        'x_BarPositions': [], #middle position of bars in x
        'x_res': [], #residual: mean-barycenter in x, per bar
        'y_res': [], #residual: mean-barycenter in y, per bar
        'QDC_perBar': [], #sum of energy loss per bar
        'lambdaX_list': [], #lambda x per bar, used to calculate lambda x of entire plane
        'lambdaY_list': [], #lambda y per bar, used to calculate lambda y of entire plane
#         'lambdaX_plane': [], 
#         'lambdaY_plane': [],
        'bary_px': [], #average momentum of barycenter in x
        'bary_pz': [], #average momentum of barycenter in z 
        "weight": [], #weights
        "relative_yBar": []
        }
        
#relative QDC method
        x_min, x_max = df_bary["coordX"].min(), df_bary["coordX"].max()
        y_min, y_max = df_bary["coordY"].min(), df_bary["coordY"].max()
        #relative method
        for b in range(amount):
            df_bary_bar = df_bary[df_bary["Bar"]==b]
            QDC = df_bary_bar["Energy_loss"].sum() #energy loss in the specified bar
#             print(QDC/QDC_tot)

            if QDC >= self.QDCbar_threshold: #energy loss threshold, usually 1MeV
                #get bar positions   
                detID_list = df_bary_bar["detectorID"].tolist() #detID of this bar
                x_right, x_left, y_top, y_bottom = barpositions[detID_list[0]]
                y_middle = (y_top+y_bottom)/2
                x_middle = (x_left+x_right)/2
                results['y_BarPositions'].append(y_middle)
                results['x_BarPositions'].append(x_middle)

                #y-barycenter:
                df_bary_x = df_bary_bar["coordX"]
                df_bary_y = df_bary_bar["coordY"]
                mean = df_bary_y.mean()
                weighted_mean = mean*(QDC/QDC_tot)
                results["relative_yBar"].append(df_bary_y*(QDC/QDC_tot))
                
                if not np.isnan(weighted_mean):
                    results['y'].append(weighted_mean)
#                     results['y_res'].append(mean - np.sum(results['y']))#does it even make sense to calculate y_res per bar?
                    
                #calculate lambda y of barycenter, if it's larger than the previous bar, then it is saved
                YMIN_local, YMAX_local = df_bary_y.min(), df_bary_y.max()
                y_max = max(y_max, df_bary_y.max())
                y_min = min(y_min, df_bary_y.min())
                lambday_bottom = (YMAX_local-YMIN_global)
                lambday_top = (YMAX_global-YMIN_local)
                results['lambdaY_list'].append(((lambday_top + lambday_bottom) * 0.5 * (QDC / QDC_tot)))

                #get average c_scint for each bar (all SiPMs), left and right
                c_scint_left, c_scint_right = [], []
                detID = str(detID_list[0])
                for i in [0, 1, 3, 4, 6, 7]:
                    detID_combined = f"{detID}_{i}"
                    if detID_combined in cscint:
                        c_scint_left.append(cscint[detID_combined][0])
                for j in [8, 9, 11, 12, 14, 15]:
                    detID_combined = f"{detID}_{j}"
                    if detID_combined in cscint:
                        c_scint_right.append(cscint[detID_combined][0])
                
                #average csint for left and right SiPMs
                c_scint_left_mean = np.mean(c_scint_left)
                c_scint_right_mean = np.mean(c_scint_right)
                c_scint_bar_mean = (c_scint_left_mean + c_scint_right_mean) / 2 #average cscint per bar
                #variance of the average cscint for left and right SiPMs
                c_scint_left_var = np.var(c_scint_left)
                c_scint_right_var = np.var(c_scint_right)

                #x-barycenter calc:
                X = df_bary.loc[(df_bary["Bar"] == b), "coordX"]
                
                #smear the signal with cscint
                t_left, t_right = X.min()/c_scint_left_mean, X.max()/c_scint_right_mean
                result = (t_right + t_left)*(c_scint_bar_mean/2)
                lambdaX_temp = (t_right - t_left)*c_scint_bar_mean

                #smear with cscint
                smear_noise_left = np.random.normal(0, np.sqrt(c_scint_left_var))
                smear_noise_right = np.random.normal(0, np.sqrt(c_scint_right_var))
                smeared_result_left = result + smear_noise_left
                smeared_result_right = result + smear_noise_right
                combined_smeared_result = (smeared_result_left + smeared_result_right) / 2

                #smear lambda with cscint
                lambdaX_smearLeft = lambdaX_temp + smear_noise_left
                lambdaX_smearRight = lambdaX_temp + smear_noise_right
                smeared_lambda = (lambdaX_smearLeft + lambdaX_smearRight)/2

                if not np.isnan(combined_smeared_result):
                    residual = X.mean() - combined_smeared_result #residual (mean - barycenter)
#                     results['x_res'].append(residual)
                    combined_smeared_result *= (QDC / QDC_tot)
                    results['x'].append(combined_smeared_result)
                    smeared_lambda *= (QDC / QDC_tot) #lambdax, smeared
                    results['lambdaX_list'].append(smeared_lambda)
                
                #get momentum of hits in barycenter
                results['bary_px'].append(df_bary.loc[df_bary["Bar"] == b, "px"].values)
                results['bary_pz'].append(df_bary.loc[df_bary["Bar"] == b, "pz"].values)
                results['QDC_perBar'].append(QDC) #QDC in a bar 
                #get weight of event
                weight_event = df_bary_bar["w"].iloc[0]
                results["weight"].append(weight_event)
                
                if df_bary_x.min()<x_min:
                    x_min = df_bary_x.min()
                if df_bary_x.max()<x_max:
                    x_max = df_bary_x.max()
                if df_bary_y.min()<y_min:
                    y_min = df_bary_y.min()
                if df_bary_y.max()<y_max:
                    y_max = df_bary_y.max()

        #the actual barycenters
        y_bary = np.sum(results['y'])
        x_bary = np.sum(results['x'])
        lambdaX = np.sum(results['lambdaX_list'])
        lambdaY = np.sum(results['lambdaY_list'])
        
        #alt lambda-y method
#         lambdaX = x_max-x_min
#         lambdaY = y_max-y_min
        
        df_bary_x, df_bary_y = df_bary["coordX"], df_bary["coordY"]
        resX = df_bary_x.mean() - x_bary
        resY = df_bary_y.mean() - y_bary
        results["x_res"].append(resX)
        results["y_res"].append(resY)
        if results["relative_yBar"]:
            flattened_yBar = np.concatenate([series.values for series in results["relative_yBar"]])
            lambdaY2 = np.sqrt(np.var(flattened_yBar))
        else:
            lambdaY2 = None
        
        #average over momentum of hits in barycenter
        if results['bary_px'] and results['bary_pz']:
            mean_px = np.mean(np.concatenate(results['bary_px']))
            mean_pz = np.mean(np.concatenate(results['bary_pz']))
        else:
            mean_px, mean_pz = 0, 0

        #return        
        if results['x'] and results['y']:
            return Result(x_bary, y_bary, results['y_BarPositions'], results['x_BarPositions'], results['x_res'], 
                          results['y_res'], results['QDC_perBar'], lambdaX, lambdaY,
                          mean_px, mean_pz, interaction_wall, results["weight"], lambdaY2)
        else:
            return None, None, None, None, None, None, None, None, None, None, None

    #residual between barycenter and mean x-&y-position of shower
    def residual(self, df_shower):
        residual_y, residual_x = [], []
        total_res_MC = []
        for e in df_shower["Eventnumber"].unique():
            for p in df_shower.loc[(df_shower["Eventnumber"] == e), "PlaneNR"].unique():
                s = 2
                empty = df_shower.loc[((df_shower["Eventnumber"] == e)&(df_shower["StationNR"] == s)
                                       &(df_shower["PlaneNR"] == p)), "coordX"].empty
                if empty == False:
                    df_temp = df_shower.loc[((df_shower["Eventnumber"] == e)&(df_shower["StationNR"] == s)
                                             &(df_shower["PlaneNR"] == p))]
                    result = self.barycenter(df_temp, p)
                    if all(entry is not None for entry in result):
                        x_bary, y_bary = result.x_bary, result.y_bary
                        if not ((y_bary == None) or (x_bary == None)):
                            X = df_temp["coordX"].mean()
                            Y = df_temp["coordY"].mean()
                            residual_y.append(Y-y_bary)
                            residual_x.append(X-x_bary)

                        elif ((y_bary == None) or (x_bary == None)):
                            pass
                    elif result is None:
                        pass
                elif empty == True:
                    pass

        #y_residuals
        #plt.figure(figsize=(10, 6))
        plt.hist(residual_y, bins=50, label="filtered y-barycenter residuals")
        plt.xlabel('y-residual (cm)')
        plt.ylabel('Counts')
        plt.title('residuals in y\n Mean - barycenter')
        plt.legend()
        plt.show()

        #x-residuals
        #plt.figure(figsize=(10, 6))
        plt.hist(residual_x, bins=50) 
        plt.xlabel('x-residual (cm)', fontsize=14)
        plt.ylabel('frequency', fontsize=14)
        plt.title('Residual = Mean - barycenter', fontsize=14)
        plt.show()
        
    #draw the bars in the xy-histogram for US
    def DrawBars(self, fig, axes):

        text_x, text_y = [], []
        USx, USy = 82.5, 6.

        if not isinstance(axes, np.ndarray):
            axes = [axes]

        # Plot each rectangle
        for detID, (Ax, Bx, Ay, By) in barpositions.items():

            plane=int(detID)%10000//1000
            # Define the coordinates of the four corners of the rectangle
            x1, y1 = Ax, Ay-USy/2  # Lower left corner
            x2, y2 = Bx, By-USy/2  # Lower right corner
            x3, y3 = Bx, By+USy/2  # Upper right corner
            x4, y4 = Ax, Ay+USy/2  # Upper left corner    
            text_x.append(x1)
            text_y.append(y1)

            for ax in axes:
                polygon_vertices = ((x1, y1), (x2, y2), (x3, y3), (x4, y4))
                polygon = Polygon(polygon_vertices, closed=True, color='grey', alpha=0.04)
                ax.add_patch(polygon)

                #ax.text((x1+x4)/2, (y1+y4)/2, str(detID%10000//1000), ha='center', va='center')


        return fig, axes

    #2d histogram in US of shower, also contains muons and reconstructed muons
    def XY_histogram_US(self, data, shower, muon, e, plot_res_perBar, process):
        df_shower_e = shower[shower["Eventnumber"] == e]
        df_muon_e = muon[muon["Eventnumber"] == e]
        p_min = 0
        p_max = 4
        Y_min, Y_max = 0, 80
        X_min, X_max = -100, 20
        x_residuals = []
        y_residuals = []

        #get clusters
#         cluster = muon_recon.make_cluster(data, e)
        df_DS = data[(data["StationNR"]==3)&(data["Eventnumber"]==e)]
        clusters = muon_recon.clustering(df_DS)

        while p_min <= p_max:
            df_shower_temp = df_shower_e.loc[df_shower_e["PlaneNR"]==p_min]
            df_muon_temp = df_muon_e.loc[df_muon_e["PlaneNR"]==p_min]
            X, Y = df_shower_temp["coordX"], df_shower_temp["coordY"]

            if not (X.empty or Y.empty):
                weights = df_shower_temp["Energy_loss"]
                #get barycenters
                result = self.barycenter(df_shower_temp, p_min, process=process)
                if all(entry is not None for entry in result):
                    lambdaX_plane, lambdaY_plane = result.lambdaX, result.lambdaY
                    x_bary, y_bary =  result.x_bary, result.y_bary
                    y_BarPositions, x_ResPerBar = result.y_BarPositions, result.x_res
                if not ((y_bary == None) or (x_bary == None)):
                    fig, ax = plt.subplots(figsize=(10, 8))

                    #plot shower
                    scatter = ax.scatter(X, Y, s=3, c=weights, cmap=plt.cm.jet, norm=mpl.colors.LogNorm())
                    #h = ax.hist2d(X, Y, bins=200, norm=mpl.colors.LogNorm(), cmap=plt.cm.jet, weights = weights)
                    cbar = fig.colorbar(scatter, ax=ax)
                    #cbar = fig.colorbar(h[3], ax=ax)
                    cbar.set_label('Energy loss (GeV)', fontsize=14)
                    plt.title(f'Event {e}\n US, Plane {p_min}', fontsize=14)
                    plt.xlabel('X coordinate (cm)', fontsize=14)
                    plt.ylabel('Y coordinate (cm)', fontsize=14)  

                    #plot residual of x-barycenter per bar for shower
                    if plot_res_perBar == True:
                        ax.axvline(0, color='black', linestyle='-', linewidth=1, alpha=0.3)
                        ax.scatter([x - min(x_ResPerBar) for x in x_ResPerBar], 
                                   [y for y in y_BarPositions if y is not None], 
                                   color="orange", label = "x-residuals per bar")

                    #plot cluster
                    label = False
                    coordX, coordY, coordZ = [], [], []
                    for cluster in clusters:
                        for hit in cluster["hits"]:
                            if label == False:
                                plt.scatter(hit["coordX"], hit["coordY"], s=50, color="red", label="reconstructed Muon " 
                                            + str(hit["pdg_code"]))
                                label = True
                            elif label == True:
                                plt.scatter(hit["coordX"], hit["coordY"], s=50, color="red")
                            X_cluster, Y_cluster, Z_cluster = hit["coordX"], hit["coordY"], hit["coordZ"]
                            coordX.append(X_cluster)
                            coordY.append(Y_cluster)
                            coordZ.append(Z_cluster)
                    points = Points([[x, y, z] for x, y, z in zip(coordX, coordY, coordZ)])
                    line_fit = Line.best_fit(points)
                    z_pos_plane = zPos["MuFilter"][20+p_min]
                    point_on_line = line_fit.point
                    dir_vector = line_fit.direction
                    t = (z_pos_plane - point_on_line[2]) / dir_vector[2]
                    muonx_plane = point_on_line[0] + t * dir_vector[0]
                    muony_plane = point_on_line[1] + t * dir_vector[1]
                    plt.scatter(muonx_plane, muony_plane, s=50, color="red", zorder=3,
                                label="reconstructed muon: (" f"{round(muonx_plane, 2)}cm, {round(muony_plane, 2)}cm)")

                    #plot muon
                    for track in df_muon_temp["TrackID"].unique():
                        X_muon = df_muon_temp.loc[(df_muon_temp["TrackID"] == track), "coordX"]
                        Y_muon = df_muon_temp.loc[(df_muon_temp["TrackID"] == track), "coordY"]
                        pdg = df_muon_temp.loc[df_muon_temp["TrackID"] == track, "pdg_code"].values[0]
                        rounded_X = round(X_muon.values[0], 2)
                        rounded_Y = round(Y_muon.values[0], 2)
#                         ax.scatter(X_muon, Y_muon, s=50, color="black", label="Muon " + str(pdg) 
#                                    + ": (" + str(rounded_X) + "cm, " + str(rounded_Y) +"cm)")

                    #labels and limits and barycenters
                    plt.axhline(y=Y.mean(), color='b', linestyle='--', linewidth=1, alpha=0.5, 
                               label=("mean: (" + str(round(X.mean(), 2)) + "cm, " + str(round(Y.mean(), 2)) + "cm)"))
                    plt.axvline(X.mean(), color='b', linestyle='--', linewidth=1, alpha=0.5)
                    xmin, xmax = x_bary-(lambdaX_plane/2), x_bary+(lambdaX_plane/2)
                    ymin, ymax = y_bary-(lambdaY_plane/2), y_bary+(lambdaY_plane/2)
                    plt.scatter(x_bary, y_bary, s=50, color="green", label=("barycenter: (" + str(round(x_bary, 2)) + "cm, "
                                                                            + str(round(y_bary, 2)) + "cm)"))
                    plt.hlines(y=y_bary, xmin=xmin, xmax=xmax, color="green", linewidth=2, label=r'$\lambda$ in x and y')
                    plt.vlines(x=x_bary, ymin=ymin, ymax=ymax, color="green", linewidth=2)
                    plt.legend(fontsize=14)
                    plt.xlim(X_min, X_max)  
                    plt.ylim(Y_min, Y_max) 

                    #plot detector bars
                    self.DrawBars(fig,ax)

                    p_min += 1

                    #append residuals to return
                    x_residuals.append(X.mean()-x_bary)  
                    y_residuals.append(Y.mean()-y_bary)
                elif ((y_bary == None) or (x_bary == None)):
                    print("Plane " + str(p_min) + ": none of the MC points pass the energy threhsold.")
                    p_min += 1
            elif (X.empty or Y.empty):
                print("Plane " + str(p_min) + ": doesn't contain any MC points.")
                p_min += 1
        plt.show() 
        #return x_residuals, y_residuals

#makes 3d tracks for muon and shower in US
class build_track:
    def __init__(self, energy_loss_threshold):
        self.threshold = energy_loss_threshold
    
    #get dimensions and positions of planes
    def planes(self, df_shower):
        p = 0
        xL, xR = [], []
        yB, yT = [], []
        while p < 5:
            df_temp = df_shower.loc[((df_shower["StationNR"] == 2)&(df_shower["PlaneNR"] == p))]
            detID1, detID2 = 20000, 20009
            x_right, x_left, y_top, y_bottom = barpositions[detID1]
            xL.append(x_left)
            xR.append(x_right)
            yB.append(y_bottom)
            x_right, x_left, y_top, y_bottom = barpositions[detID2]
            yT.append(y_top)
            p += 1

        #z positions of the US planes, first digit is subsystem, second digit is planeNR
        zPos
        plane = 20
        zPos_muon = []
        while plane < 25:
            zPos_muon.append(zPos["MuFilter"][plane])
            plane += 1
        plane = 30
        while plane < 37:
            if plane % 2 == 0:
                zPos_muon.append(zPos["MuFilter"][plane])
            plane += 1
           
        #make xL, xR, yT, yB the same length as zPos_muon by adding None values
        max_len = max(len(xL), len(xR), len(yT), len(yB), len(zPos_muon))
        xL.extend([None] * (max_len - len(xL)))
        xR.extend([None] * (max_len - len(xR)))
        yT.extend([None] * (max_len - len(yT)))
        yB.extend([None] * (max_len - len(yB)))

        plane_pos = list(zip(xL, xR, yT, yB, zPos_muon))
        return plane_pos
    
    def reconstruct_muon(self, df):
        df_DS = df[df["StationNR"]==3]
        muon_recon = muon_reconstruction()
        clusters = muon_recon.clustering(df_DS)
        coordX, coordY, coordZ = [], [], []
        muonX, muonY, muonZ = [], [], []
        for cluster in clusters:
            for hit in cluster["hits"]:
                X_cluster, Y_cluster, Z_cluster = hit["coordX"], hit["coordY"], hit["coordZ"]
                coordX.append(X_cluster)
                coordY.append(Y_cluster)
                coordZ.append(Z_cluster)
        if (coordX and coordY and coordZ):
            coordinates = list(zip(coordX, coordY, coordZ))
            coordinates.sort(key=lambda coord: coord[2])  # Sort by the Z coordinate
            coordX_sorted, coordY_sorted, coordZ_sorted = zip(*coordinates)
            points = Points([[x, y, z] for x, y, z in zip(coordX_sorted, coordY_sorted, coordZ_sorted)])
#             points = Points([[x, y, z] for x, y, z in zip(coordX, coordY, coordZ)])
            line_fit = Line.best_fit(points)
        else:
            points, line_fit = None, None
            
        return points, line_fit, clusters
        
    #make barycenter for muon&shower, all events
    #take average z-position per plane
    def make_points(self, data, df_shower, df_muon, process):

        df_DS = data[data["StationNR"]==3]
        df_DS = df_DS[df_DS["pz"]>100]
        clusters = muon_recon.clustering(df_DS)
        coordX, coordY, coordZ = [], [], []
        muonX, muonY, muonZ = [], [], []
        for cluster in clusters:
            for hit in cluster["hits"]:
                X_cluster, Y_cluster, Z_cluster = hit["coordX"], hit["coordY"], hit["coordZ"]
                coordX.append(X_cluster)
                coordY.append(Y_cluster)
                coordZ.append(Z_cluster)
        if (coordX and coordY and coordZ):
            points = Points([[x, y, z] for x, y, z in zip(coordX, coordY, coordZ)])
            line_fit = Line.best_fit(points)
            for p in range(0, 5):
                z_pos_plane = zPos["MuFilter"][20+p]
                point_on_line = line_fit.point
                dir_vector = line_fit.direction
                t = (z_pos_plane - point_on_line[2]) / dir_vector[2]
                muonx_plane = point_on_line[0] + t * dir_vector[0]
                muony_plane = point_on_line[1] + t * dir_vector[1]
                muonX.append(muonx_plane)
                muonY.append(muony_plane)
                muonZ.append(z_pos_plane)
        else:
            muonX, muonY, muonZ = None, None, None            

        #create shower barycenters for each plane, all events
        #mean z position per plane
        showerX_bary, showerY_bary, showerZ = [], [], []
        showerY_QDC, BarsFired, y_Bars = [], [], []
        df_shower_US = df_shower.loc[df_shower["StationNR"]==2]
        for p in df_shower_US["PlaneNR"].unique():
            result = bary_setup.barycenter(df_shower, p, process="numu")
            if all(entry is not None for entry in result):
                showerX_bary.append(result.x_bary)
                showerY_bary.append(result.y_bary)
                showerY_QDC.append(result.QDC_perBar)
                y_Bars.append(result.y_BarPositions)
                z = df_shower_US.loc[(df_shower_US["PlaneNR"]==p), "coordZ"].mean()
                if not np.isnan(z):
                    showerZ.append(z)

        return muonX, muonY, muonZ, showerX_bary, showerY_bary, showerZ, showerY_QDC, y_Bars, clusters
    
    #plot shower and muon in each US plane, 2d
    def plot_perPlane(self, e, muon_points, shower_points, y_Bars, QDC, lambdaX_plane, lambdaY_plane, baryX, 
                      baryY, planeNR):
        y_Bars = [y for y in y_Bars if y is not None]
        length = len(y_Bars)
        length = len(QDC)
        showerX_list = [shower_points[0]] * length
        fig, ax = plt.subplots(figsize=(6, 6))
        plt.scatter(muon_points[0], muon_points[1], s=50, c='red', label="reconstructed muon")
        sc = plt.scatter(showerX_list, y_Bars, s=50, c=QDC, label="shower barycenters")
        #plot clusters
        label = False
#         for cluster in clusters:
#             for hit in cluster["hits"]:
#                 if label == False:
#                     plt.scatter(hit["coordX"], hit["coordY"], color="orange", label="Cluster " + str(hit["TrackID"]) 
#                                + "\n(pdg: " + str(hit["pdg_code"]) + ")")
#                     label = True
#                 elif label == True:
#                     plt.scatter(hit["coordX"], hit["coordY"], color="orange")
        
        #plot lambda of barycenter
        xmin, xmax = baryX-(lambdaX_plane/2), baryX+(lambdaX_plane/2)
        ymin, ymax = baryY-(lambdaY_plane/2), baryY+(lambdaY_plane/2)
        plt.hlines(y=baryY, xmin=xmin, xmax=xmax, color="green", linewidth=1, label=r'$\lambda$ in x and y')
        plt.vlines(x=baryX, ymin=ymin, ymax=ymax, color="green", linewidth=1)
        
        plt.title("Event " + str(e) + ", Plane " + str(planeNR), fontsize=14)
        plt.xlabel('X (cm)', fontsize=14)
        plt.ylabel('Y (cm)', fontsize=14)
        plt.xlim(-100, 20)
        plt.ylim(0, 80)
        plt.legend(fontsize=14)
        cbar = fig.colorbar(sc, ax=ax)
        cbar.set_label(label='Energy loss (GeV)\n sum over bar', fontsize=14)
        bary_setup.DrawBars(fig,ax)
        plt.tight_layout()
        plt.show()
        
    #get angle between shower and muon tracks
    def get_angle(self, line1, line2):
        direction1 = np.array(line1.direction)
        direction2 = np.array(line2.direction)
        dot_product = np.dot(direction1, direction2)
        angle_radians = np.arccos(dot_product)
        angle_degrees = np.degrees(angle_radians)
        if dot_product < 0:
            angle_degrees = 180 - angle_degrees

        #return angle_degrees
        return angle_radians
    
    #make and plot tracks for shower and muon in US, 3d
    def track_3d(self, e, df, df_shower, df_muon, plot_3d, plot_planes, polar, azimuth, process):        
        muonX, muonY, muonZ, showerX, showerY, showerZ, QDC_perBar, y_Bars, clusters = track.make_points(df,
        df_shower, df_muon, process=process)
        if not all([muonX, muonY, muonZ, showerX, showerY, showerZ]):
            print("Points couldn't be made")
            return None, None, None, None, None

        #Muon
        points_muon = Points([[x, y, z] for x, y, z in zip(muonX, muonY, muonZ)])
        #check for collinear points
#         if points_muon.are_collinear(tol=0.0001):
#             print("Muon points are collinear")
#             return None, None, None, None, None
        line_fit_muon = Line.best_fit(points_muon)
        direction_muon = line_fit_muon.direction.unit()
        direction_components_muon = direction_muon.to_array()
        polar_angle_muon = np.arccos(direction_components_muon[2]) #polar angle of muon
        azimuthal_angle_muon = np.arctan2(direction_components_muon[1], direction_components_muon[0]) #azimuthal angle of muon
        #polar_angle = np.arccos(direction.z)

        #shower
        points_shower = Points([[x, y, z] for x, y, z in zip(showerX, showerY, showerZ)])
        #check for collinear points
#         if points_shower.are_collinear(tol=0.0001):
#             print("Shower points are collinear")
#             return None, None, None, None, None
        line_fit_shower = Line.best_fit(points_shower)
        angle = self.get_angle(line_fit_muon, line_fit_shower) #angle between muon and shower track in radians
        direction_shower = line_fit_shower.direction.unit()
        direction_components_shower = direction_shower.to_array()
        polar_angle_shower = np.arccos(direction_components_shower[2]) #polar angle of shower
        azimuthal_angle_shower = np.arctan2(direction_components_shower[1], direction_components_shower[0]) #azimuthal angle of shower
        #polar_angle2 = np.arccos(direction2.z)
        
        amount_planes = min(points_muon.shape[0], points_shower.shape[0])
        
        def plot_square(ax, x_left, x_right, y_top, y_bottom, z):
            vertices = np.array([
                [x_left, y_top, z],
                [x_right, y_top, z],
                [x_right, y_bottom, z],
                [x_left, y_bottom, z],
                [x_left, y_top, z]  
            ])

            square = Poly3DCollection([vertices], alpha=0.08, linewidth=1, facecolor='black')
            ax.add_collection3d(square)

        def plot_label(ax, x, y, z, p):
            text = "Plane " + str(p)
            ax.text(x, y, z, text, color='black')
        
        if plot_3d == True:
            fig = plt.figure(figsize=(10, 8))
            ax = fig.add_subplot(111, projection='3d')
            fig.subplots_adjust(left=0.1)

            points_muon.plot_3d(ax, c='red', s=70, depthshade=False)
            line_fit_muon.plot_3d(ax, t_1=-50, t_2=50, c='red', label="reconstructed muon track", alpha=0.8)
            points_shower.plot_3d(ax, c='blue', s=70, depthshade=False, label="shower barycenters")
            line_fit_shower.plot_3d(ax, t_1=-50, t_2=50, c='blue', label="shower track", alpha=0.8)
            
                
            #get cluster data
            X_cluster, Y_cluster, Z_cluster = [], [], []
            counter = 1
            for cluster in clusters:
                for hit in cluster["hits"]:
#                     ax.scatter(hit["coordX"], hit["coordY"], hit["coordZ"], label="Cluster " + str(counter) 
#                                + "\npdg: " + str(hit["pdg_code"]))
                    counter += 1  
                    X_cluster.append(hit["coordX"])
                    Y_cluster.append(hit["coordY"])
                    Z_cluster.append(hit["coordZ"])

            if (X_cluster and Y_cluster and Z_cluster):
                points_cluster = Points([[x, y, z] for x, y, z in zip(X_cluster, Y_cluster, Z_cluster)])
                line_fit_cluster = Line.best_fit(points_cluster) 
                points_cluster.plot_3d(ax, c='green', s=70, depthshade=False, label="cluster per plane")
                line_fit_cluster.plot_3d(ax, t_1=-50, t_2=50, c='green', label="reconstructed muon track", alpha=0.8)

            ax.set_xlim([-100, 20])
            ax.set_ylim([0, 80])
            #ax.set_zlim([370, 475])
            ax.set_zlim([370, 550])
            ax.set_title("Event " + str(e), fontsize=14)
            ax.set_xlabel('X', fontsize=14)
            ax.set_ylabel('Y', fontsize=14)
            ax.set_zlabel('Z', fontsize=14)
            ax.view_init(elev=polar, azim=azimuth)
            #ax.legend()
            ax.legend(loc='center left', fontsize=12, bbox_to_anchor=(-0.3, 0.15))
            half_pi = np.pi / 2
            if angle > half_pi:
                angle = np.pi - angle
            if polar_angle_muon > half_pi:
                polar_angle_muon = np.pi - polar_angle_muon
            if polar_angle_shower > half_pi:
                polar_angle_shower = np.pi - polar_angle_shower
            ax.text2D(0.05, 0.95, "Angles:\n" + r"$\alpha$ = " + f"{round(np.degrees(angle), 2)}°\n" 
                      + r"$\theta_{shower}$ = " + f"{round(np.degrees(polar_angle_shower), 2)}°\n" + r"$\theta_{muon} = $" 
                      + f"{round(np.degrees(polar_angle_muon), 2)}°", transform=ax.transAxes, 
                      fontsize=12, verticalalignment='top')

            plane_pos = self.planes(df)

            p = 0
            while p < 5:
                plot_square(ax, plane_pos[p][0], plane_pos[p][1], plane_pos[p][2], plane_pos[p][3], plane_pos[p][4])
                plot_label(ax, plane_pos[p][0], plane_pos[p][2], plane_pos[p][4], p)
                p += 1
            p = 5
            while p < len(plane_pos):
                plot_square(ax, plane_pos[0][0], plane_pos[0][1], plane_pos[0][2], plane_pos[0][3], plane_pos[p][4])
                plot_label(ax, plane_pos[0][0], plane_pos[0][2], plane_pos[p][4], p)
                p += 1

            plt.show()

        if plot_planes == True:
            p = 0
            while p < amount_planes:
                result = bary_setup.barycenter(df_shower, p, process=process)
                lambdaX_plane, lambdaY_plane = result.lambdaX, result.lambdaY
                baryX, baryY = result.x_bary, result.y_bary
                self.plot_perPlane(e, points_muon[p,:2], points_shower[p,:2], y_Bars[p], QDC_perBar[p], lambdaX_plane, lambdaY_plane,
                                   baryX, baryY, planeNR = p)
                p += 1

        return angle, polar_angle_muon, azimuthal_angle_muon, polar_angle_shower, azimuthal_angle_shower
    
#muon reconstruction and cuts on muon momentum and slope
class muon_reconstruction:
    def __init(self):
        pass
    
    def clustering(self, data):
        # Filter the data
#         data = data[data["pz"] > 130]
#         if data.empty:
#             print("what")
        data = data.reset_index(drop=True)
        #print(data)
        clusters = []
        hitDict = {}

        # Populate hitDict with valid hits
        for index, row in data.iterrows():
            detectorID = row['detectorID']
            if (detectorID // 10000) < 3:
                continue
            hitDict[detectorID] = index

        # Create hit list
        hitList = list(hitDict.keys())
        if len(hitList) > 0:
            hitList.sort()  # Smallest detID first: sorted low z to high z

            tmp = [hitList[0]]  # Temporary list to append c (detID) if c is a neighbor of cprev
            cprev = hitList[0]  # c previous
            ncl = 0  # Number of clusters in plane
            last = len(hitList) - 1

            for i in range(1, len(hitList)):  # Start from second element
                c = hitList[i]
                neighbour = False
                if (c - cprev) == 1 or (c - cprev) == 2:
                    neighbour = True  # Allow for one missing channel to be a neighbor
                    tmp.append(c)
                if not neighbour or c == hitList[last] or c % 1000 == 59:  # Cluster is done for this hit if one of the conditions is fulfilled. Last condition corresponds to last bar
                    first = tmp[0]  # First element of cluster
                    N = len(tmp)  # Amount of hits in cluster
                    hitvector = []

                    for aHit in tmp:
#                         print(len(hitDict) , aHit)
                        hitvector.append(data.iloc[hitDict[aHit]])

                    # Define cluster object, add to cluster list
                    aCluster = {
                        "detID": first,
                        "N": N,
                        "hits": hitvector,
                        "neighbour": neighbour
                    }
                    clusters.append(aCluster)

                    if c != hitList[last]:
                        ncl += 1  # Add to amount of clusters
                        tmp = [c]  # Starts next cluster
                    elif not neighbour:  # Save last channel
                        hitvector = []
                        hitvector.append(data.iloc[hitDict[c]])

                        aCluster = {
                            "detID": c,
                            "N": 1,
                            "hits": hitvector,
                            "neighbour": neighbour
                        }
                        clusters.append(aCluster)
                cprev = c

        return clusters
    
    def reconstruct_muon(self, df):
        df_DS = df[df["StationNR"]==3]
        clusters = self.clustering(df_DS)
        coordX, coordY, coordZ = [], [], []
        muonX, muonY, muonZ = [], [], []
        polar_angle = None
        if len(clusters) > 2:
            for cluster in clusters:
                for hit in cluster["hits"]:
                    X_cluster, Y_cluster, Z_cluster = hit["coordX"], hit["coordY"], hit["coordZ"]
                    coordX.append(X_cluster)
                    coordY.append(Y_cluster)
                    coordZ.append(Z_cluster)

            if coordX and coordY and coordZ:
                # Sort coordinates by Z to ensure correct order
                coordinates = list(zip(coordX, coordY, coordZ))
                coordinates.sort(key=lambda coord: coord[2])  # Sort by the Z coordinate
                coordX_sorted, coordY_sorted, coordZ_sorted = zip(*coordinates)

                # Perform linear fits for X and Y against Z
                slope_x, intercept_x, _, _, _ = linregress(coordZ_sorted, coordX_sorted)
                slope_y, intercept_y, _, _, _ = linregress(coordZ_sorted, coordY_sorted)

                # Define a function for the fitted line
                line_fit_x = lambda z: slope_x * z + intercept_x
                line_fit_y = lambda z: slope_y * z + intercept_y
                
                polar_angle = math.atan(math.sqrt(slope_x**2 + slope_y**2))

                # Package the fits in a dictionary for convenience
                line_fit = {
                    "slope_x": slope_x, "intercept_x": intercept_x,
                    "slope_y": slope_y, "intercept_y": intercept_y,
                    "line_fit_x": line_fit_x, "line_fit_y": line_fit_y
                }

                points = coordinates

            else:
                points, line_fit, polar_angle = None, None, None
        else: points, line_fit, polar_angle = None, None, None
            
        return points, line_fit, clusters, polar_angle
    
    def reconstruction_residual(self, df_muon, df, delta_value):
        diff_dict_x = {0: [], 1: [], 2: [], 3: []}
        diff_dict_y = {0: [], 1: [], 2: [], 3: []}
        chi_squaredX_list, chi_squaredY_list = [], []
        df_muon = df_muon[df_muon["TrackID"]==1]
        
        for e in df_muon["Eventnumber"].unique():
            df_muon_e = df_muon[df_muon["Eventnumber"]==e]
            df_e = df[df["Eventnumber"]==e]

            for p in range(0, 4):
                df_muon_e_plane = df_muon_e[df_muon_e["PlaneNR"]==p]
                muonX, muonY = df_muon_e_plane["coordX"].tolist(), df_muon_e_plane["coordY"].tolist()
                
                df_e_plane = df_e[df_e["PlaneNR"]==p]
                
                plane_offset = p if p == 0 else p + 1
                DS_plane = zPos["MuFilter"][30 + plane_offset]
                    
                muon_track = muon_recon.reconstruct_muon(df_e_plane) 
                if muon_track is None:
                    continue
                track_points = muon_track[0]
#                 track_vector = muon_track[1]
#                 if track_points is not None:
#                     for point in track_points:
#                         x_0, y_0, z_0 = point
#                         v_x, v_y, v_z = track_vector.direction
#                         t = (DS_plane - z_0) / v_z
#                         trackX = x_0 + t * v_x
#                         trackY = y_0 + t * v_y
#                         diff_listX = [trackX - muon for muon in muonX]
#                         diff_listY = [trackY - muon for muon in muonY]
#                         for diff_x in diff_listX:
#                             diff_dict_x[p].append(diff_x)
#                         for diff_y in diff_listY:
#                             diff_dict_y[p].append(diff_y)

                line_fit = muon_track[1]
                if track_points is not None and line_fit is not None:
                    trackX = line_fit["line_fit_x"](DS_plane)
                    trackY = line_fit["line_fit_y"](DS_plane)

                    # Calculate the difference between reconstructed track and muon hit coordinates
                    diff_listX = [trackX - muon for muon in muonX]
                    diff_listY = [trackY - muon for muon in muonY]

                    # Append differences to corresponding dictionaries for each plane
                    for diff_x in diff_listX:
                        diff_dict_x[p].append(diff_x)
                    for diff_y in diff_listY:
                        diff_dict_y[p].append(diff_y)


        # Now calculate standard deviations (uncertainties) per plane just once
#         uncertainty_X = {p: np.std(diff_dict_x[p]) for p in range(4)}
#         uncertainty_Y = {p: np.std(diff_dict_y[p]) for p in range(4)}
        
        delta = delta_value  # or any appropriate value for your use case
        uniform_std = delta / np.sqrt(3)  # standard deviation for a uniform distribution from -delta to +delta

        # Set uniform uncertainties per plane
        uncertainty_X = {p: uniform_std for p in range(4)}
        uncertainty_Y = {p: uniform_std for p in range(4)}

        total_dof_X = 0
        total_dof_Y = 0
        # Calculate chi-squared in one pass after knowing the uncertainties
        for p in range(0, 4):
            for diff_x in diff_dict_x[p]:
                chi_squaredX = (diff_x ** 2) / uncertainty_X[p] ** 2
                chi_squaredX_list.append(chi_squaredX)

            for diff_y in diff_dict_y[p]:
                chi_squaredY = (diff_y ** 2) / uncertainty_Y[p] ** 2
                chi_squaredY_list.append(chi_squaredY)
        
            num_points_X = len(diff_dict_x[p])
            num_points_Y = len(diff_dict_y[p])
            dof_X = num_points_X - 2
            dof_Y = num_points_Y - 2
            total_dof_X += dof_X
            total_dof_Y += dof_Y
            
        reduced_chi_squared_X = sum(chi_squaredX_list) / total_dof_X if total_dof_X > 0 else None
        reduced_chi_squared_Y = sum(chi_squaredY_list) / total_dof_Y if total_dof_Y > 0 else None

        return diff_dict_x, diff_dict_y, chi_squaredX_list, chi_squaredY_list, reduced_chi_squared_X, reduced_chi_squared_Y
                    
    def hit_distribution(self, df):
        hitBars = []
        for e in df["Eventnumber"].unique():
            df_e = df.loc[df["Eventnumber"]==e]
            clusters = self.clustering(df_e)
            
            for cluster in clusters:
                for hit in cluster["hits"]:
                    plane = hit["PlaneNR"]
                    bar  = hit["Bar"]
                    hitBars.append([plane, bar])
                    
        planes = {0: [], 1: [], 2: [], 3: []}

        for entry in hitBars:
            plane_index = entry[0]
            if plane_index in planes:
                planes[plane_index].append(entry[1])

        fig, axs = plt.subplots(4, sharex=True, figsize=(8, 6))
        default_colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
        bins = np.arange(0, 120, 1)
        
        for key in planes:
            axs[key].hist(planes[key], bins=bins, color=default_colors[key], label=f"Plane {key}")
            axs[key].set_ylabel("Counts", fontsize=14)
            axs[key].tick_params(axis='x',labelbottom=True)
            axs[key].legend(fontsize=14)
        
        axs[3].set_xlabel("Bar", fontsize=14)
        axs[0].set_title("distribution of hit bars from clustering algorithm", fontsize=14)
        plt.tight_layout()  # Adjust layout for better spacing
        plt.show()
        
        return hitBars

class FiducialCutProcessor:
    def __init__(self, vertical_min_ds, vertical_max_ds, horizontal_min_ds, horizontal_max_ds,
                 vertical_min_scifi, vertical_max_scifi, horizontal_min_scifi, horizontal_max_scifi,
                 reversed_ds, reversed_scifi):
        self.vertical_min_ds = vertical_min_ds
        self.vertical_max_ds = vertical_max_ds
        self.horizontal_min_ds = horizontal_min_ds
        self.horizontal_max_ds = horizontal_max_ds
        self.reversed_ds = reversed_ds
        
        self.vertical_min_scifi = vertical_min_scifi
        self.vertical_max_scifi = vertical_max_scifi
        self.horizontal_min_scifi = horizontal_min_scifi
        self.horizontal_max_scifi = horizontal_max_scifi
        self.reversed_scifi = reversed_scifi

    def avg_ds_fiducial_cut(self, df):
        avg_ver = 0.
        n_ver = 0
        avg_hor = 0.
        n_hor = 0

        for _, hit in df.iterrows():
            if hit['StationNR'] == 3:
                x = hit['detectorID'] % 1000
                isvertical = hit["Bar"] > 59

                if isvertical:
                    avg_ver += x
                    n_ver += 1
                else:
                    avg_hor += x
                    n_hor += 1

        if (n_ver + n_hor) == 0:
            plot_var = [-1, -1]
            return False

        avg_ver = avg_ver / n_ver if n_ver else -1
        avg_hor = avg_hor / n_hor if n_hor else -1

        if n_ver == 0 or n_hor == 0:
            return False

        if avg_hor < self.horizontal_min_ds or avg_hor > self.horizontal_max_ds or \
           avg_ver < self.vertical_min_ds or avg_ver > self.vertical_max_ds:
            return False

        return True

    def avg_scifi_fiducial_cut(self, df):
        avg_ver = 0.
        n_ver = 0
        avg_hor = 0.
        n_hor = 0

        for _, hit in df.iterrows():
            mat = int(hit["detectorID"] / 10000) % 10
            sipm = int(hit["detectorID"] / 1000) % 10
            channel = hit["detectorID"] / 1000
            x = channel + sipm * 128 + mat * 4 * 128

            isvertical = int(hit["detectorID"] / 100000) % 10 == 1

            if isvertical:
                avg_ver += x
                n_ver += 1
            else:
                avg_hor += x
                n_hor += 1

        if (n_ver + n_hor) == 0:
            plot_var = [-1, -1]
            return False

        avg_ver = avg_ver / n_ver if n_ver else -1
        avg_hor = avg_hor / n_hor if n_hor else -1

        if n_ver == 0 or n_hor == 0:
            return False

        if not self.reversed_scifi:
            if avg_hor < self.horizontal_min_scifi or avg_hor > self.horizontal_max_scifi or \
               avg_ver < self.vertical_min_scifi or avg_ver > self.vertical_max_scifi:
                return False
        else:
            if self.horizontal_min_scifi < avg_hor < self.horizontal_max_scifi or \
               self.vertical_min_scifi < avg_ver < self.vertical_max_scifi:
                return False

        return True
    
    def apply_cuts(self, df_hcal, df_scifi, plot):
        event_numbers = df_hcal["Eventnumber"].unique()

        # Initialize counters
        amount_pass_hcal = 0
        amount_pass_scifi = 0

        # Keep track of events to be kept
        valid_events_hcal = set()
        valid_events_scifi = set()

        for event_number in event_numbers:
            event_data_hcal = df_hcal[df_hcal["Eventnumber"] == event_number]
            event_data_scifi = df_scifi[df_scifi["Eventnumber"] == event_number]

            pass_scifi = self.avg_scifi_fiducial_cut(event_data_scifi)
            pass_hcal = self.avg_ds_fiducial_cut(event_data_hcal)

            if pass_hcal:
                amount_pass_hcal += 1
                valid_events_hcal.add(event_number)
            if pass_scifi:
                amount_pass_scifi += 1
                valid_events_scifi.add(event_number)

        # Combine valid event sets to get all valid events
        valid_events_combined = valid_events_hcal.intersection(valid_events_scifi)
        
        # Filter DataFrames based on valid events
        df_hcal_filtered = df_hcal[df_hcal["Eventnumber"].isin(valid_events_combined)]
        df_scifi_filtered = df_scifi[df_scifi["Eventnumber"].isin(valid_events_combined)]

        # For plotting purposes, compute the filtered DataFrames for specific conditions
        df_cutScifi = df_hcal[df_hcal["Eventnumber"].isin(valid_events_scifi)]
        df_cutDS = df_hcal[df_hcal["Eventnumber"].isin(valid_events_hcal)]

        if plot:
            self.plot_cuts(df_hcal, df_cutDS, df_cutScifi, df_hcal_filtered)
            
        return df_hcal_filtered, df_scifi_filtered, amount_pass_hcal, amount_pass_scifi

    def plot_cuts(self, df, df_cutDS, df_cutScifi, df_filtered):    
        # Define possible values for the bins
        possible_values = np.arange(60, 121)

        # Data for pre-cut histogram
        pre_cut = df.loc[((df["StationNR"] == 3) & (df["Bar"] > 59)), "Bar"]

        # Create histograms from the raw data
        fig = plt.figure(figsize=(12, 6))
        ax = fig.add_subplot(111)

        # Plot pre-cut histogram
        ax.hist(pre_cut, bins=possible_values, histtype="step", linewidth=1.5, label='Pre-cut')

        # Filtered data
        filtered_data_DS = df_cutDS.loc[((df_cutDS["StationNR"] == 3) & (df_cutDS["Bar"] > 59)), "Bar"]
        filtered_data_Scifi = df_cutScifi.loc[((df_cutScifi["StationNR"] == 3) & (df_cutScifi["Bar"] > 59)), "Bar"]
        post_cut = df_filtered.loc[((df_filtered["StationNR"] == 3) & (df_filtered["Bar"] > 59)), "Bar"]

        # Plot filtered histograms
        ax.hist(filtered_data_Scifi, bins=possible_values, histtype="step", linewidth=1.5, label='Scifi cut')
        ax.hist(filtered_data_DS, bins=possible_values, histtype="step", linewidth=1.5, label='DS cut')
#         ax.hist(post_cut, bins=possible_values, histtype="bar", alpha=0.4, label='Post fiducial cuts')

        # Labels and title
        ax.set_xlabel('BarNR', fontsize=14)
        ax.set_ylabel('Frequency', fontsize=14)
        ax.set_title('Frequency of hits in DS vertical bars\n after avg DS&SciFi fiducial cuts', fontsize=14)
        ax.legend(fontsize=14) 

        # Add cut criteria to the plot
        cut_name = "Avg DS Ver bar in [70,105], \nHor in [10,50]; \nAvg SciFi Ver channel in [200,1200], \nHor in [300,1336]"
        textbox_content = "Cut criteria:\n" + cut_name
        ax.text(130, 0, textbox_content, fontsize=14, bbox=dict(facecolor='white', alpha=0.5, pad=3))
        
        plt.tight_layout()
        plt.show()

class Tracking():
    def Init(self):
        self.systemAndPlanes = {1:2, 2:5, 3:7} #amount of planes per station
        self.DSnPlanes = 3 #amount of planes required for a track
        self.DSnHits = 5 #restrict max nr. of clusters per plane
        millimeter = 1./10.
        self.sigmaMufiUS_spatial = 2.*10.*millimeter #resolution US
        self.sigmaMufiDS_spatial = 0.3*10.*millimeter #resolution DS
        self.Debug = False
        self.fitter = ROOT.genfit.KalmanFitter()
        self.fitter.setMaxIterations(50)
        fM = ROOT.genfit.FieldManager.getInstance()
        bfield = ROOT.genfit.ConstField(0,0,0)   # constant field of zero
        fM.init(bfield)
        geoMat = ROOT.genfit.TGeoMaterialInterface()
        ROOT.genfit.MaterialEffects.getInstance().init(geoMat)
        ROOT.genfit.MaterialEffects.getInstance().setNoEffects()

        # if self.simulation: self.snd_geo = SndlhcGeo.GeoInterface(options.path + options.geoFile)
        # else: 
        #     if path.find('eos')>0: path  = options.server+options.path         
        #     self.snd_geo = SndlhcGeo.GeoInterface(path+options.geoFile)

        filepath = "/eos/experiment/sndlhc/MonteCarlo/MuonBackground/muons_up/scoring_2.5/geofile_full.Ntuple-TGeant4.root"
        self.snd_geo = SndlhcGeo.GeoInterface(filepath)
        self.MuFilter = self.snd_geo.modules['MuFilter']

        # self.Scifi    = self.snd_geo.modules['Scifi']
        # self.systemAndPlanes = {1:2,2:5,3:7}
        # self.zPos = self.getAverageZpositions()

    def clustering(self, data):
        data = data.reset_index(drop=True)
        clusters = []
        hitDict = {}

        # Populate hitDict with valid hits
        for index, row in data.iterrows():
            detectorID = row['detectorID']
            if (detectorID // 10000) < 3:
                continue
            hitDict[detectorID] = index

        # Create hit list
        hitList = list(hitDict.keys())
        if len(hitList) > 0:
            hitList.sort()  # Smallest detID first: sorted low z to high z

            tmp = [hitList[0]]  # Temporary list to append c (detID) if c is a neighbor of cprev
            cprev = hitList[0]  # c previous
            ncl = 0  # Number of clusters in plane
            last = len(hitList) - 1

            for i in range(1, len(hitList)):  # Start from second element
                c = hitList[i]
                neighbour = False
                if (c - cprev) == 1 or (c - cprev) == 2:
                    neighbour = True  # Allow for one missing channel to be a neighbor
                    tmp.append(c)
                if not neighbour or c == hitList[last] or c % 1000 == 59:  # Cluster is done for this hit if one of the conditions is fulfilled. Last condition corresponds to last bar
                    first = tmp[0]  # First element of cluster
                    N = len(tmp)  # Amount of hits in cluster
                    hitvector = []

                    for aHit in tmp:
                        hitvector.append(data.iloc[hitDict[aHit]])

                    # Define cluster object, add to cluster list
                    aCluster = {
                        "detID": first,
                        "N": N,
                        "hits": hitvector,
                        "neighbour": neighbour
                    }
                    clusters.append(aCluster)

                    if c != hitList[last]:
                        ncl += 1  # Add to amount of clusters
                        tmp = [c]  # Starts next cluster
                    elif not neighbour:  # Save last channel
                        hitvector = []
                        hitvector.append(data.iloc[hitDict[c]])

                        aCluster = {
                            "detID": c,
                            "N": 1,
                            "hits": hitvector,
                            "neighbour": neighbour
                        }
                        clusters.append(aCluster)
                cprev = c
        return clusters
        
    def ExecuteTask(self, df, nPlanes=3):
        self.trackCandidates = {}
        self.clusMufi.Delete()
        self.dsCluster()
#         self.trackCandidates['DS'] = self.DStrack()
        i_muon = -1
        
        hcal_grouped = df.groupby("Eventnumber")
        for event, group in hcal_grouped:
            muon_points, muon_track, clusters, polar_angle = muon_recon.reconstruct_muon(group)
            trackcandidates = tracking.DStrack(event, clusters)
            for candidate in trackcandidates:
                for key, value in candidate.items():
                    rc = tracking.fitTrack(value)
                    
                    if type(rc)==type(1):
                        # rc==-2: not converged, rc==-1 not consistent
                        print('trackfit failed',rc,aTrack)
                    else:
                        print('track fitted apparently')
                        i_muon += 1
                        rc.SetUniqueID(3)

                        # add the tracks
                        if self.genfitTrack:
                            self.fittedTracks.Add(rc)
                        else:
                            # Load items into snd track class object
                            #if not rc.getFitStatus().isFitConverged(): continue
                            this_track = ROOT.sndRecoTrack(rc)
                            pointTimes = ROOT.std.vector(ROOT.std.vector('float'))()
                            for pnt in rc.getPointsWithMeasurement():
                                hitID = pnt.getRawMeasurement().getHitId()
                                cluster = self.clusMufi[hitID]
                                pointTimes.push_back([cluster.GetTime()])
                            this_track.setRawMeasTimes(pointTimes)
                            this_track.setTrackType(rc.GetUniqueID())

                            # Store the track in sndRecoTrack format
                            self.fittedTracks[i_muon] = this_track
        
    def DStrack(self, event, clusters):
        trackCandidates = []

        stations = {}
        clusPerStation = {}
        planesPerProjection = {0: 0, 1: 0}

        s = 3
        for p in range(self.systemAndPlanes[s]+1): 
            stations[s*10+p] = {}
            clusPerStation[s*10+p] = 0
            
        k=-1
        for cluster in clusters:
            k+=1
            detID = cluster["detID"]
            if detID//10000 < 3:
                continue # only work with DS hits
            p = (detID//1000)%10 
            bar = detID%1000
            plane = s*10+p
            if bar < 60: 
                plane = s*10+2*p  # even numbers horizontal planes, odd numbers vertical planes
            else:  
                plane = s*10+2*p+1
            stations[plane][k] = cluster #assign cluster to corresponding plane in dict
            clusPerStation[plane] += 1 #increment cluster count per station
            for p in clusPerStation:
                if clusPerStation[p]>0:
                    planesPerProjection[p%2]+=1 # Count number of horizontal (p%2==0) and vertical (p%2==1) hits

        failed = False
        #require 3 planes to be hit for a track, otherwise return empty list
        if planesPerProjection[1]<self.DSnPlanes or planesPerProjection[0]<self.DSnPlanes: 
            return trackCandidates

        # require max 1 cluster per plane if only 2 planes hit per projection, otherwise return empty list
        if self.DSnPlanes==2 and (planesPerProjection[1]==2 or planesPerProjection[0]==2):            
            for p in clusPerStation:
                if clusPerStation[p]>1:
                    failed = True
                    break
            if failed: 
                return trackCandidates
            hitlist = {}
            for p in stations:
                for k in stations[p]:  
                    hitlist[k] = stations[p][k]
                trackCandidates.append(hitlist)
                return trackCandidates

        #if more than 5 clusters per plane, return empty list
        for p in clusPerStation:
            if clusPerStation[p]>self.DSnHits: 
                return trackCandidates

        # proj = 0, horizontal, max 3 planes
        # proj = 1, vertex,     max 4 planes
        # require a plane with 1 cluster as seed
        seed = -1
        combinations = {}
        hitlist = {}
        #seed plane selection
        for proj in range(2): #horizontal and vertical projection
            for plane in range(self.systemAndPlanes[s]+1): #loop through all planes in detector
                if not plane%2==proj: 
                    continue
                if clusPerStation[s*10+plane]==1: #suitable plane is found
                    seed = s*10+plane # seed defined as the lowest z-position plane to have 1 cluster! 
                    break
            if seed < 0: #if no suitable seed is found then exit function
                return trackCandidates

            #track building 
            A, B = ROOT.TVector3(), ROOT.TVector3()
            combinations[proj] = {}  #empty dictionary for storing combinations of hits per projection
            for keyA in stations[seed]: #get average cluster position of seed plane
                clusterA = stations[seed][keyA]
                detID = clusterA["detID"]
                self.MuFilter.GetPosition(detID, A, B)
                posA = (A+B)/2.
                
                #search for clusters in other planes
                clInPlane = {}
                for p2 in range(self.systemAndPlanes[s]+1): #for each projection
                    if not p2%2==proj: 
                        continue
                    planeB = s*10+p2
                    if planeB == seed: #skip plane of the seed
                        continue
                    l = len(stations[planeB])
                    if l>0: #if plane contains clusters
                        clInPlane[planeB] = l
                srt = sorted(clInPlane) #sort planes

                #check total number of hits in the two planes closest to seed plane
                if len(srt) < 2:
                    return trackCandidates
                if len(stations[srt[0]])+len(stations[srt[1]]) > self.DSnHits+1: 
                    return trackCandidates
                planeB = srt[0] #plane closest to seed plane for analysis
                for keyB in stations[planeB]:
                    clusterB = stations[planeB][keyB]
                    detID = clusterB["detID"]
                    self.MuFilter.GetPosition(detID, A, B)
                    posB = (A+B)/2.

                    delBA = posB-posA #displacement vector from the seed cluster to the cluster in plane B
                    if proj==0: #x = slope*z+b
                        slope = delBA[1]/delBA[2] #vertical displacement/z
                        b = posA[1]-slope*posA[2] #intercept of track in plane
                    else: #y = slope*z+b
                        slope = delBA[0]/delBA[2] #horizontal displacement/z
                        b = posA[0]-slope*posA[2]

                    #next plane for analysis
                    for p3 in range(self.systemAndPlanes[s]+1):
                        if not p3%2==proj: 
                            continue
                        planeC = s*10+p3
                        if planeC == seed or planeC == planeB: #not seed plane ot plane B
                            continue
                        for keyC in stations[planeC]:
                            clusterC = stations[planeC][keyC]
                            detID = clusterC["detID"]
                            self.MuFilter.GetPosition(detID, A, B)
                            posC = (A+B)/2.

                            eX = posC[2]*slope+b #expected position
                            if proj==0: 
                                res = abs(eX-posC[1])
                            else:       
                                res = abs(eX-posC[0])
                            if proj==0: 
                                combinations[proj][res] = [[seed,keyA], [planeB,keyB], [planeC,keyC]]
                            else:
                                #fourth plane
                                for p4 in range(self.systemAndPlanes[s]+1):
                                    if not p4%2==proj: 
                                        continue
                                    planeD = s*10+p4
                                    if planeD == seed or planeD == planeB or planeD == planeC: 
                                        continue
                                    if len(stations[planeD])==0:
                                        combinations[proj][res] = [[seed,keyA], [planeB,keyB], [planeC,keyC]]
                                    for keyD in stations[planeD]:
                                        clusterD = stations[planeD][keyD]
                                        detID = clusterD["detID"]
                                        self.MuFilter.GetPosition(detID, A, B)
                                        posD = (A+B)/2.

                                        eX = posD[2]*slope+b
                                        if proj==0: 
                                            res+= abs(eX-posD[1])
                                        else:       
                                            res+= abs(eX-posD[0])
                                        combinations[proj][res] = [[seed, keyA], [planeB, keyB], [planeC, keyC], [planeD, keyD]]

            # find combination with smallest residual
            srt = sorted(combinations[proj])[0]
            for x in combinations[proj][srt]:
                hitlist[x[1]] = stations[x[0]][x[1]]

        trackCandidates.append(hitlist)
        return trackCandidates
    
    def fitTrack(self, event, hitlist):
        # hitlist:  clusterID: [A,B] endpoints of scifiCluster
        hitPosLists = {}
        trID = 0

        posM = ROOT.TVector3(0, 0, 0.)
        momM = ROOT.TVector3(0, 0, 100.)  # default track with high momentum
        
        # approximate covariance
        covM = ROOT.TMatrixDSym(6)
        res = self.sigmaMufiDS_spatial #resolution of DS 
        for  i in range(3):   covM[i][i] = res*res 
        for  i in range(3,6): covM[i][i] = ROOT.TMath.Power(res / (4.*2.) / ROOT.TMath.Sqrt(3), 2)
        # Runge-Kutta track representation(q/p, u', v', u, v). # Runge-Kutta track representation(q/p, u', v', u, v).
        # u,v: positions on a genfit.DetPlane. 13 for muon.
        rep = ROOT.genfit.RKTrackRep(13) 
        
        # start state
        state = ROOT.genfit.MeasuredStateOnPlane(rep)
        rep.setPosMomCov(state, posM, momM, covM)
        
        # create track
        seedState = ROOT.TVectorD(6)
        seedCov   = ROOT.TMatrixDSym(6)
        rep.get6DStateCov(state, seedState, seedCov)
        # initialize a track object using a seed state (initial estimate of parameters)
        theTrack = ROOT.genfit.Track(rep, seedState, seedCov)
        
        # make measurements sorted in z
        unSortedList = {}
        tmpList = {}
        A, B = ROOT.TVector3(), ROOT.TVector3() #create 2 empty vectors
        i = 0
        for cluster in hitlist:
            # print(f"Processing cluster: {event, i}")
            # print(cluster)
            detID = cluster["detID"]
            detSys = 3    

            self.MuFilter.GetPosition(detID, A, B)

            distance = 0
            tmp = array('d', [A[0], A[1], A[2], B[0], B[1], B[2], distance])
            unSortedList[A[2]] = [ROOT.TVectorD(7,tmp), detID, cluster, detSys]

        # print(unSortedList)
        sorted_z=list(unSortedList.keys())
        sorted_z.sort()
        if len(sorted_z) < 3:
            return -1
        for z in sorted_z:
            tp = ROOT.genfit.TrackPoint() 
            hitCov = ROOT.TMatrixDSym(7)
            res = self.sigmaMufiUS_spatial
            maxDis = 5.0
            hitCov[6][6] = res*res

            measurement = ROOT.genfit.WireMeasurement(unSortedList[z][0], hitCov, 1, 6, tp) # the measurement is told which trackpoint it belongs to
            measurement.setMaxDistance(maxDis)
            measurement.setDetId(unSortedList[z][1])
            # measurement.setHitId(unSortedList[z][2])
            tp.addRawMeasurement(measurement) # package measurement in the TrackPoint                                  
            theTrack.insertPoint(tp)  # add point to Track
        if not theTrack.checkConsistency():
            print("track not consistent")
            theTrack.Delete()
            return -1

        if self.Debug:
            print(f"Track initialized with {theTrack.getNumPoints()} points.")
            for p in theTrack.getPointsWithMeasurement():
                rawM = p.getRawMeasurement()
                print(f"Point DetID: {rawM.getDetId()}, Raw Hit Coords: {rawM.getRawHitCoords()}")

        
        # do the fit
        self.fitter.processTrack(theTrack) # resortHits bool=False by default.
        # processTrack(theTrack) is inherited from genfit::AbsFitter. 
        # Uses genfit.processTrackWithRep(track, trackRep, false) with all hits
        # processTrack(theTrack) starts with the cardinalRep (lowest chi2?) then proceeds to others without sorting.
        # How is cardinal rep chosen? Does it matter?
        # From Thomas: processTrackWithRep(theTrack,rep,True)... from Andrew: bool resortHits cannot be passed as True? 
        # Andrew: according to https://eic.github.io/doxygen/d0/dba/classgenfit_1_1AbsFitter.html#a03e279c67ca889f5fe6480eb8a1691cf)

        fitStatus = theTrack.getFitStatus()

        if not fitStatus.isFitConverged():
            # print("Track fit did not converge.")
            return -1
        if self.Debug: 
            print("Fit result: converged chi2 Ndf", fitStatus.isFitConverged(), fitStatus.getChi2(), fitStatus.getNdf())
        # if not fitStatus.isFitConverged() and 0 > 1:
        #     theTrack.Delete()
        #     return -1
        if self.Debug: 
            chi2 = fitStatus.getChi2()/(fitStatus.getNdf()+1E-15)
            fittedState = theTrack.getFittedState()
            P = fittedState.getMomMag()
            print("track fitted Ndf #Meas P",fitStatus.getNdf(), theTrack.getNumPointsWithMeasurement(),P)
            for p in theTrack.getPointsWithMeasurement():
                rawM = p.getRawMeasurement()
                info = p.getFitterInfo()
                if not info: 
                    continue
                detID = rawM.getDetId()
                print(detID,"weights",info.getWeights()[0],info.getWeights()[1],fitStatus.getNdf())

        return theTrack