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



def extract_us_signal(ch, N):
    signal_sum = 0
    us_points = {"Eventnumber": [], "px": [], "py": [], "pz": [], "detectorID": [], "time": [], "pdg_code": [], "Energy": [], "Energy_loss": [], "coordX": [], "coordY": [], "coordZ": []}
    for hit in ch.MuFilterPoint:   

        station = int(hit.GetDetectorID()/1000000)
        # if station == 0:
        #     print(station)
        P = np.sqrt(hit.GetPx()**2 + hit.GetPy()**2 + hit.GetPz()**2)
        E = np.sqrt(hit.GetPx()**2 + hit.GetPy()**2 + hit.GetPz()**2 + ch.MCTrack[hit.GetTrackID()].GetMass()**2)
        pdg = hit.PdgCode()
        if pdg in [22, 111, 113, 2112]:
            continue
        us_points["detectorID"].append(hit.GetDetectorID())
        time = hit.GetTime()
        # if time > 25. or time < 0.:
        #     continue
        signal_sum += hit.GetEnergyLoss()
        us_points["px"].append(hit.GetPx())
        us_points["py"].append(hit.GetPy())
        us_points["pz"].append(hit.GetPz())
        us_points["time"].append(hit.GetTime())
        us_points["pdg_code"].append(hit.PdgCode())
        us_points["Energy"].append(E)
        us_points["Energy_loss"].append(hit.GetEnergyLoss())
        us_points["coordX"].append(hit.GetX())
        us_points["coordY"].append(hit.GetY())
        us_points["coordZ"].append(hit.GetZ())
        us_points["Eventnumber"].append(N)
    return us_points

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
    parser.add_argument(
        "--genie",
        dest="genie",
        type=bool,
        default=False)

    args = parser.parse_args()
    f = r.TFile.Open(args.outputfile, 'recreate')
    h = {}
    f.cd()
    ch = r.TChain('cbmsim')
    eos = "root://eosuser.cern.ch/"
    energy = args.Energy
    genie = args.genie
    if not genie:
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
    else:
        basePath = ["/eos/experiment/sndlhc/MonteCarlo/Neutrinos/Genie/sndlhc_13TeV_down_volMuFilter_20fb-1_SNDG18_02a_01_000/1/sndLHC.Genie-TGeant4.root"]
    for base in basePath:
        # print(base)
        ch.Add(str(base))
    # Define histograms

    # hist_list = {l: r.TH1I("plane" + f"_{l}", "plane" + f"_{l}; {label};", 200, bin_min, bin_max)}
    #ch.Add(args.inputfile)
    print(ch.GetListOfBranches())
    data_general = {}
    for N in range(args.nStart, args.nStart + args.nEvents):
        ch.GetEvent(N)
        if N % 1000 == 0:
            print(f"Event {N}")
            us_signal = extract_us_signal(ch, N)
            if not len(data_general):
                for key in us_signal:
                    data_general[key] = us_signal[key]
                    #print(key, len(data_general[key]))

            for key in us_signal:
                data_general[key] += us_signal[key]

            #scifi_signal = extract_scifi_signal(ch)
            #print(us_signal)
        
    # print(wall_info)
    # vis_info(wall_info)
    # vis_info_parts(wall_info)

    #print(signal_data)
    data_general = pd.DataFrame(data_general)
    data_general.to_csv(f"/eos/user/t/tismith/SWAN_projects/genie_ana_output/us_output.csv")
    #print(signal_data.shape)
    #signal_relation(signal_data, energy)
    #reco_resol(signal_data, energy)
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
