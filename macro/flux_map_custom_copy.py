#!/usr/bin/env python2
from __future__ import division
import argparse
import numpy as np
import ROOT as r
# Fix https://root-forum.cern.ch/t/pyroot-hijacks-help/15207 :
r.PyConfig.IgnoreCommandLineOptions = True
import shipunit as u
import rootUtils as ut
import logger as log
from array import array
from pathlib import Path
import pandas as pd
from math import floor
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns 
sns.set_style("whitegrid")

PDG_conv = {11: "e-", -11: "e+", 2212: "p", 211: "pi+", -211: "pi-", 1000060120: "C", 321: "K+", -321: "K-", 1000020040: "Ca", 13: "mu-", -13: "mu+"}

def vis_info(wall_info):
    fig, ax = plt.subplots(3, 4, figsize = (16, 12), dpi = 150)
    for i, wall in enumerate(wall_info):
        part_nums = np.array([[len(event[st]) for st in event] for event in wall_info[wall]]).T
        for j, part_hist_info in enumerate(part_nums):
            ax[i][j].grid()
            ax[i][j].hist(part_hist_info, bins = 100, label = "$\\overline{N_{p}} = " + "{}$".format(int(part_hist_info.mean())), color = "red")
            ax[i][j].set_title(f"Shower starts in wall {i+1}. Station {j+1}")
            ax[i][j].legend()
            
    fig.savefig("part_output.pdf")
   
 
def vis_info_parts(wall_info):
    fig, ax = plt.subplots(3, 4, figsize = (16, 12), dpi = 150)
    part_list = {i:{j: {} for j in range(1, 5)} for i in range(1,4)}
    for i, wall in enumerate(wall_info):
        for event in wall_info[wall]:
            for j, st in enumerate(event):
                unique, counts = np.unique(event[st], return_counts=True)
                for freq in zip(list(unique), list(counts)):
                    if freq[0] not in part_list[i+1][j+1]:
                        part_list[i+1][j+1][freq[0]] = [freq[1]]
                    else:
                        part_list[i+1][j+1][freq[0]].append(freq[1])

    for i in range(3):
        for j in range(4): 
            part_list_sorted_av = {key: np.mean(part_list[i+1][j+1][key]) for key in part_list[i+1][j+1]}
            part_list_sorted_av_sorted = dict(sorted(part_list_sorted_av.items(), key=lambda x:x[1], reverse=True)[:5])
            if part_list_sorted_av == {}:
                continue
            # print(part_list_sorted_av_sorted)        
            ax[i][j].grid()
            ax[i][j].bar(list(map(lambda x: PDG_conv[x], list(part_list_sorted_av_sorted.keys()))), [part_list_sorted_av_sorted[key] for key in part_list_sorted_av_sorted], color = "red")
            ax[i][j].set_title(f"Shower starts in wall {i+1}. Station {j}")
            # ax[i][j].legend()
            
    fig.savefig("part_output_parts.pdf")
    
def isVertical(detid):
    if floor(detid/100000)%10 == 1: return True
    else: return False
def read_files(chain, eos, filepath, filename, file_num):
    for i in range(1, file_num + 1):
        chain.Add(eos + filepath + "/" + str(i) + "/" + filename)



def extract_us_signal(ch):
    signal_sum = 0
    for hit in ch.MuFilterPoint:   
        station = int(hit.GetDetectorID()/1000000)
        # if station == 0:
        #     print(station)
        P = np.sqrt(hit.GetPx()**2 + hit.GetPy()**2 + hit.GetPz()**2)
        pdg = hit.PdgCode()
        if pdg in [22, 111, 113, 2112]:
            continue
        time = hit.GetTime()
        if time > 25. or time < 0.:
            continue
        signal_sum += hit.GetEnergyLoss() 
    return signal_sum
def extract_scifi_signal(ch):
    signal_sum = 0
    for hit in ch.ScifiPoint:   
        station = int(hit.GetDetectorID()/1000000)
        if station == 0:
            print(station)
        P = np.sqrt(hit.GetPx()**2 + hit.GetPy()**2 + hit.GetPz()**2)
        pdg = hit.PdgCode()
        if pdg in [22, 111, 113, 2112]:
            continue
        time = hit.GetTime()
        if time > 25 or time < 0:
            continue
        signal_sum += hit.GetEnergyLoss() 
    return signal_sum        

def signal_relation(signal_data_init, energy):
    fig, ax = plt.subplots(2, 2, figsize = (8,8), dpi = 200)

    h = ax[1][1].hist2d(signal_data_init["signal_us"], signal_data_init["signal_scifi"], norm=mpl.colors.LogNorm(), bins = (50,50),
                         cmap = plt.cm.jet)
    fig.colorbar(h[3], ax = ax[1][1])    
    wall = 1
    for i in range(2):
        for j in range(2):
            if i == 1 and j == 1:
                ax[i][j].set_title("All")
                continue
            signal_data = signal_data_init.query(f"wall == {wall}").copy()
            h = ax[i][j].hist2d(signal_data["signal_us"], signal_data["signal_scifi"], norm=mpl.colors.LogNorm(), bins = (50,50),
                            cmap = plt.cm.jet)
            fig.colorbar(h[3], ax = ax[i][j])
            ax[i][j].set_title(f"Shower starts at wall {wall}")
            wall += 1
    ax[0][0].set_ylabel("Scifi signal sum/event [GeV]")
    ax[1][0].set_ylabel("Scifi signal sum/event [GeV]")
    ax[1][0].set_xlabel("US signal sum/event [GeV]")
    ax[1][1].set_xlabel("US signal sum/event [GeV]")
    fig.savefig(f"signal_distribution_{energy}.pdf")
    
    
def reco_resol(Data, energy):
    fig, ax = plt.subplots(figsize = (8,8), dpi = 200)
    ax.hist(Data["reco"], bins = 100)
    ax.set_xlabel("Energy [GeV]")
    ax.set_title(f"Pion energy {energy} GeV")
    fig.savefig(f"reco_{energy}.pdf")
def main():
    parser = argparse.ArgumentParser()
    #parser.add_argument(
    #    'inputfile',
    #    help='''Simulation results to use as input. '''
    #    '''Supports retrieving files from EOS via the XRootD protocol.''')
    parser.add_argument(
        '-o',
        '--outputfile',
        default='flux_map.root',
        help='''File to write the flux maps to. '''
        '''Will be recreated if it already exists.''')
    # parser.add_argument(
    # '-P',
    # '--Pcut',
    # default= 0.0,
    # help='''set momentum cut''')
    # parser.add_argument(
    # '-E',
    # '--Eloss',
    # default= 0.0,
    # help= '''set Eloss cut''')
    parser.add_argument(
        '--nStart',
        dest="nStart",
        type = int,
        default=0)
    
    parser.add_argument(
        '--nEvents',
        dest="nEvents",
        type = int,
        default=1)
    
    parser.add_argument(
        '--Energy',
        dest="Energy",
        type = int,
        default=100)  

    parser.add_argument(
        '--Merged',
        dest="merged",
        type = bool,
        default=True)     
    args = parser.parse_args()
    f = r.TFile.Open(args.outputfile, 'recreate')
    h = {}
    f.cd()
    ch = r.TChain('cbmsim')
    eos = "root://eosuser.cern.ch/"
    energy = args.Energy
    if not args.merged:
        if energy != 300:
            filepath = f"/eos/user/e/ekhaliko/Documents/SND_Data/test_{energy}GeV_n10k_aug2023_pi+/"
            fileName = "sndLHC.PG_211-TGeant4.root"
        else:
            filepath = f"/eos/user/e/ekhaliko/Documents/SND_Data/test_{energy}GeV_n10k_aug2023_pi-/"
            fileName = "sndLHC.PG_-211-TGeant4.root"
        # filename = "sndLHC.PG_-211-TGeant4.root"
        path = filepath
        basePath = sorted(Path(path).glob(f'**/{fileName}'))
        print("{} files to read in {}".format(len(basePath), path))
    else:
        if energy != 300:
            filepath = f"/eos/user/e/ekhaliko/Documents/SND_Data/test_{energy}GeV_n10k_aug2023_pi+/"
            fileName = "merge.root"
        else:
            filepath = f"/eos/user/e/ekhaliko/Documents/SND_Data/test_{energy}GeV_n10k_aug2023_pi-/"
            filepath = "/eos/user/e/ekhaliko/Documents/SND_Data/test_300GeV_n100k_aug2023_pi-_new"
            fileName = "merge.root"        
        path = filepath
        basePath = sorted(Path(path).glob(f'{fileName}'))
        print("{} files to read in {}".format(len(basePath), path))
    
    for base in basePath:
        # print(base)
        ch.Add(str(base))
    # Define histograms

    # hist_list = {l: r.TH1I("plane" + f"_{l}", "plane" + f"_{l}; {label};", 200, bin_min, bin_max)}
    #ch.Add(args.inputfile)
    print(ch.GetListOfBranches())
    N_plane_ZY = {i: 0 for i in range(1, 5)}
    N_plane_ZX = {i: 0 for i in range(1, 5)}
    N_plane_ZY_part = {i: [] for i in range(1, 5)}
    N_plane_ZX_part = {i: [] for i in range(1, 5)}
    wall_info = {i:[] for i in range(1,4)}
    event_freq = {i:0 for i in range(5)}
    signal_data = {"event_number": [], "signal_us": [], "signal_scifi": [], "reco": [], "wall": []}
    for N in range(args.nStart, args.nStart + args.nEvents):
        ch.GetEvent(N)
        if N % 1000 == 0:
            print(f"Event {N}")
        mc_hits = 0
        #print(P)
        # check_1wall, check_2wall, check_3wall = [],[],[]
        # true_1wall, true_2wall, true_3wall = 1, 1, 1
        N_plane_ZY = {i: 0 for i in range(1, 5)}
        N_plane_ZX = {i: 0 for i in range(1, 5)}
        N_plane_ZX_mom = {i: 0 for i in range(1, 5)}
        N_plane_ZY_part = {i: [] for i in range(1, 5)}
        N_plane_ZX_part = {i: [] for i in range(1, 5)}
        track_list = {i: {} for i in range(1, 5)} 
        track_list_mom = {i: {} for i in range(1, 5)} 
        station_start = 0
        k_start = 0
        for hit in ch.ScifiPoint:
            station = int(hit.GetDetectorID()/1000000)
            if station == 0:
                print(station)
            if hit.GetTrackID() < 0:
                continue
            pdg = hit.PdgCode()
            if pdg in [22, 111, 113, 2112]:
                continue
            # print(hit.GetTrackID(), pdg, hit.GetPx(), hit.GetPy(), hit.GetPz())
            E = np.sqrt(hit.GetPx()**2 + hit.GetPy()**2 + hit.GetPz()**2 + ch.MCTrack[hit.GetTrackID()].GetMass()**2)
            # time = hit.GetTime()
            # if time > 25 or time < 0:
                # continue
            #print(hit.GetEnergyLoss())
            if isVertical(hit.GetDetectorID()):
                if hit.GetTrackID() not in track_list[station]:
                    track_list[station][hit.GetTrackID()] = 1
                    track_list[station][hit.GetTrackID()] = E
                    #print(N, station, hit.PdgCode(), hit.GetTrackID(), ch.MCTrack[hit.GetTrackID()].GetMotherId(), ch.MCTrack[hit.GetTrackID()].GetProcName(), P)
                else:
                    continue
                 
                # if hit.PdgCode() == -211 and hit.GetTrackID() == 0:
                #     if P < 299.8 and not k_start:
                #         station_start = station
                #         k_start = 1
                        
                # print(N, station, hit.PdgCode(), hit.GetTrackID(), P) 
                # if station == 2:
                #     check_1wall += []
                #     if len(check_1wall) > 1:
                #         true_1wall = 0
                # if station == 3:
                #     check_2wall += []
                # if station == 4:
                #     check_3wall += 1
                N_plane_ZX[station] += 1
                # if hit.PdgCode() not in N_plane_ZX_part[station]:
                #     N_plane_ZX_part[station][hit.PdgCode()] = []
                N_plane_ZX_mom[station] += E
                N_plane_ZX_part[station].append(hit.PdgCode())
            else:
                N_plane_ZY[station] += 1
                # if hit.PdgCode() not in N_plane_ZY_part[station]:
                #     N_plane_ZY_part[station][hit.PdgCode()] = []
                # P = np.sqrt(hit.GetPx()**2 + hit.GetPy()**2 + hit.GetPz()**2)
                N_plane_ZY_part[station].append(hit.PdgCode())
        
        # print(N_plane_ZX, N_plane_ZY)
        if not station_start:
            for station in range(1,4):
                #N_plane_ZX[station+1] > 3*N_plane_ZX[station] and 
                if 0.75*N_plane_ZX_mom[station] > N_plane_ZX_mom[station+1]:
                    # if -211 in N_plane_ZX_part[station]:
                    station_start = station
                    break
        # if station_start != 0:
        #     station_start -= 1                
        if station_start != 0 and station_start < 4:
            # print(N_plane_ZX_part)
            wall_info[station_start].append(N_plane_ZX_part) 
            signal_data["event_number"].append(N)
            us_signal = extract_us_signal(ch)
            scifi_signal = extract_scifi_signal(ch)
            signal_data["signal_us"].append(us_signal)
            signal_data["signal_scifi"].append(scifi_signal)
            signal_data["reco"].append(145*us_signal + 500*scifi_signal)        
            signal_data["wall"].append(station_start)             
        event_freq[station_start if station_start != 4 else 0] += 1
        
    # print(wall_info)
    # vis_info(wall_info)
    # vis_info_parts(wall_info)
    print(event_freq)
    #print(signal_data)
    signal_data = pd.DataFrame(signal_data)
    signal_data.to_csv(f"output_signal_{energy}_no-time_wo-neg-trid_F-alg_en-loss.csv")
    #print(signal_data.shape)
    signal_relation(signal_data, energy)
    reco_resol(signal_data, energy)
    # Data_X = pd.DataFrame(N_plane_ZX_part)
    # Data_Y = pd.DataFrame(N_plane_ZY_part)
    # print(Data_X, Data_Y)
    # Data_X.to_csv("output_X.csv")
    # Data_Y.to_csv("output_Y.csv")


    f.Close()
    #np.savetxt('test_out', np.array(Eloss_total), fmt='%f') 
    #f = open(args.norm + "/flux", "w")
    #f.write(str(B_ids_unw) + "\t" + str(B_ids))
    #f.close()


if __name__ == '__main__':
    r.gErrorIgnoreLevel = r.kWarning
    r.gROOT.SetBatch(True)
    main()
