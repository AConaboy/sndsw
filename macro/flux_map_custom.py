#!/usr/bin/env python2
from __future__ import division
import argparse,os,sys
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
import SndlhcGeo
import subprocess
sns.set_style("whitegrid")

PDG_conv = {11: "e-", -11: "e+", 2212: "p", 211: "pi+", -211: "pi-", 1000060120: "C", 321: "K+", -321: "K-", 1000020040: "Ca", 13: "mu-", -13: "mu+"}

parser = argparse.ArgumentParser()
# parser.add_argument('inputfile',help='Simulation results to use as input. Supports retrieving files from EOS via the XRootD protocol.', required=False)
parser.add_argument('-o','--outputfile',default='flux_map.root',help='File to write the flux maps to. Will be recreated if it already exists.')
parser.add_argument('-P','--Pcut',default= 0.0,help='set momentum cut')
parser.add_argument('-E','--Eloss',default= 0.0,help= 'set Eloss cut')
parser.add_argument('--nStart',dest="nStart",type = int,default=0)
parser.add_argument('--nEvents', dest="nEvents",type = int,default=10)
parser.add_argument('--Energy',dest="Energy",type = int,default=100)  
parser.add_argument('--Merged',dest="merged",type = bool,default=True)  
parser.add_argument("--genie",dest="genie",type=bool,default=False)
parser.add_argument('--inputfile',dest="inputfile",type=str,default=None,help='GENIE input file for convert_script.C')
parser.add_argument('-C', "--HTCondor", dest="HTCondor",action='store_true')
# parser.add_argument("--muonDIS", type=str, default=None)
parser.add_argument('--muonDIS', action='store_true', help='flag to indicate muon DIS')
parser.add_argument('--muonPassing', action='store_true', help='flag to indicate passing muons')
# parser.add_argument("--muonPassing", type=str, default=None)
parser.add_argument('--neutralhadrons', action='store_true', help='flag to indicate neutral hadrons')
# parser.add_argument("--neutralhadrons", type=str, default=None)

# Set eos outpath
eos_paths = {'tismith':'/eos/user/t/tismith/SWAN_projects/genie_ana_output/' , 'aconsnd':'/eos/user/a/aconsnd/SWAN_projects/Simulation/data/'} # can add Eduard
found=False
for key in eos_paths: 
    if os.getcwd().find(key)!=-1:
        eospath=eos_paths[key]
        found=True
        break
if not found:
    print('who is the user?')
    sys.exit()

args = parser.parse_args()

def GetData(event, N):
    us_signal = extract_us_signal(ch, N)
    scifi_signal = extract_scifi_signal(ch, N)
    return us_signal, scifi_signal

def EventLoop():
    us_data_general = {}
    scifi_data_general = {}

    if args.nEvents==-1:
        start, end = 0, ch.GetEntries() 
    else: start, end = args.nStart, args.nStart + args.nEvents
    # for N in range(args.nStart, args.nStart + args.nEvents):
    for N in range(start, end):
        
        #ch.GetEvent(N)
        event = ch.GetEntry(N)
        if N % 100 == 0:
            print(f"Event {N}")
        
        # Get data for us and scifi
        us_signal, scifi_signal = GetData(event, N)

        # Store US data for event N in overall dictionary
        if not len(us_data_general):
            for key in us_signal:
                us_data_general[key] = us_signal[key]
        for key in us_signal:
            us_data_general[key] += us_signal[key]

        # Store Scifi data for event N in overall dictionary
        if not len(scifi_data_general):
            for key in scifi_signal:
                scifi_data_general[key] = scifi_signal[key]
        for key in scifi_signal:
            scifi_data_general[key] += scifi_signal[key]

    return us_data_general, scifi_data_general             

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

def run_convert_script(input_file, output_file):
    cmd = f'root -l -q \'convert_script.C("{input_file}", "{output_file}")\\'
    subprocess.call(cmd, shell=True)

def read_convert_script_output(output_file):
    data = {}
    max_length = 0
    with open(output_file, 'r', encoding='latin-1') as f:
        for line in f:
            line = line.strip()

            if not line or not line[0].isdigit():
                    continue

            tokens = line.split()

            event_number = int(tokens[0])
            hit_id = int(tokens[1])
            pdg = int(tokens[2])
            first_mother = int(tokens[3])
            name = tokens[4]
            px = float(tokens[5])
            py = float(tokens[6])
            pz = float(tokens[7])
            energy = float(tokens[8])
            is_stable = bool(int(tokens[9]))
            vx = float(tokens[10])
            vy = float(tokens[11])
            vz = float(tokens[12])

            # Create a dictionary entry for each event number if it doesn't exist
            if event_number not in data:
                data[event_number] = []

            # Append the particle information to the event entry
            data[event_number].append({
                'hit_id': hit_id,
                'pdg': pdg,
                'first_mother': first_mother,
                'name': name,
                'px': px,
                'py': py,
                'pz': pz,
                'energy': energy,
                'is_stable': is_stable,
                'vx': vx,
                'vy': vy,
                'vz': vz
            })

            # Track the maximum length encountered
            current_length = len(data[event_number])
            if current_length > max_length:
                max_length = current_length

        # Fill in missing entries with None
        for event_number in data:
            while len(data[event_number]) < max_length:
                data[event_number].append({
                    'hit_id': None,
                    'pdg': None,
                    'first_mother': None,
                    'name': None,
                    'px': None,
                    'py': None,
                    'pz': None,
                    'energy': None,
                    'is_stable': None,
                    'vx': None,
                    'vy': None,
                    'vz': None
                })
    return data

def GetInteractionWall(x, y, z):
    nav = r.gGeoManager.GetCurrentNavigator()
    node = nav.FindNode(x, y, z)
    volume = node.GetVolume()

    return volume.GetName()

def extract_us_signal(ch, N):
    us_points = {"Wall": [], "start_pos" : [], "MotherID" : [], "Process": [], "TrackID" : [], "track_px": [], 
    "track_py": [], "track_pz": [], "Eventnumber": [], "px": [], "py": [], "pz": [], "detectorID": [], "time": [], 
    "pdg_code": [], "Energy_loss": [], "coordX": [], "coordY": [], "coordZ": [], "w": [], "rapidity": []}

    for hit in ch.MuFilterPoint:   
        us_points["detectorID"].append(hit.GetDetectorID())
        us_points["TrackID"].append(hit.GetTrackID())
        us_points["px"].append(hit.GetPx())
        us_points["py"].append(hit.GetPy())
        us_points["pz"].append(hit.GetPz())
        us_points["time"].append(hit.GetTime())
        us_points["pdg_code"].append(hit.PdgCode())
        us_points["Energy_loss"].append(hit.GetEnergyLoss())
        us_points["coordX"].append(hit.GetX())
        us_points["coordY"].append(hit.GetY())
        us_points["coordZ"].append(hit.GetZ())
        us_points["Eventnumber"].append(N)

        trackID = hit.GetTrackID()
        if trackID >= 0:
            mctrack = ch.MCTrack[trackID]
            mother_pdg = ch.MCTrack[ch.MCTrack[hit.GetTrackID()].GetMotherId()].GetPdgCode()
            us_points["track_px"].append(mctrack.GetPx())
            us_points["track_py"].append(mctrack.GetPy())
            us_points["track_pz"].append(mctrack.GetPz())
            us_points["MotherID"].append(mother_pdg)
            us_points["Process"].append(mctrack.GetProcName())
            us_points["start_pos"].append(mctrack.GetStartZ())
            us_points["w"].append(mctrack.GetPx())
            us_points["rapidity"].append(mctrack.GetPx())
                    
        elif trackID < 0:
            us_points["track_px"].append(0)
            us_points["track_py"].append(0)
            us_points["track_pz"].append(0)
            us_points["MotherID"].append(0)
            us_points["Process"].append(0)
            us_points["start_pos"].append(0)

    nav = r.gGeoManager.GetCurrentNavigator()
    for x, y, z in zip(us_points["coordX"], us_points["coordY"], us_points["coordZ"]):
        node = nav.FindNode(x, y, z)
        volume = node.GetVolume()
        us_points["Wall"].append(volume.GetName())

    return us_points

def extract_scifi_signal(ch, N):
    signal_sum = 0
    scifi_points = {"Eventnumber": [], "time": [], "coordX": [], "coordY": [], "coordZ": [], "TrackID": [], "px": [], "py": [], "pz": [], "detectorID": [], "pdg_code": [], "Energy_loss": []}
    for hit in ch.ScifiPoint:   
        station = int(hit.GetDetectorID()/1000000)
        if station == 0:
            print(station)

        signal_sum += hit.GetEnergyLoss() 
        scifi_points["detectorID"].append(hit.GetDetectorID())
        scifi_points["TrackID"].append(hit.GetTrackID())
        scifi_points["pdg_code"].append(hit.PdgCode())
        scifi_points["Eventnumber"].append(N)
        scifi_points["Energy_loss"].append(hit.GetEnergyLoss())
        scifi_points["px"].append(hit.GetPx())
        scifi_points["py"].append(hit.GetPy())
        scifi_points["pz"].append(hit.GetPz())
        scifi_points["time"].append(hit.GetTime())
        scifi_points["coordX"].append(hit.GetX())
        scifi_points["coordY"].append(hit.GetY())
        scifi_points["coordZ"].append(hit.GetZ())

    return scifi_points      

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

def getMuonData(ntuple):

    Muon_points = {"run": [], "event": [], "id": [], "generation": [], "E": [], "w": [], "x": [], "y": [], "px": [], "py": [], "t": [], "z": [], "pz": []}

    branches = ntuple.GetListOfBranches()
    for branch in branches:
        branch_name = branch.GetName()
        if branch_name != "W":
            Muon_points[branch_name] = []

            for entry in range(ntuple.GetEntries()):
                ntuple.GetEntry(entry)
                value = ntuple.GetLeaf(branch_name).GetValue()
                Muon_points[branch_name].append(value)

    df_Muon_points = pd.DataFrame(Muon_points)

    return df_Muon_points

def get_background_data(ch):
    data_points = {"Eventnumber": [], "px": [], "py": [], "pz": [], "detectorID": [], "time": [], "pdg_code": [], "Energy_loss": [], "coordX": [], "coordY": [], "coordZ": []}

    start, end = 0, ch.GetEntries() 
    for N in range(start, end):
        event = ch.GetEntry(N)
        if N % 1000 == 0:
            print(f"Event {N}")

    #for event in ch:
        for hit in ch.MuFilterPoint:    
            data_points["Eventnumber"].append(N)
            data_points["detectorID"].append(hit.GetDetectorID())
            data_points["px"].append(hit.GetPx())
            data_points["py"].append(hit.GetPy())
            data_points["pz"].append(hit.GetPz())
            data_points["time"].append(hit.GetTime())
            data_points["pdg_code"].append(hit.PdgCode())
            data_points["Energy_loss"].append(hit.GetEnergyLoss())
            data_points["coordX"].append(hit.GetX())
            data_points["coordY"].append(hit.GetY())
            data_points["coordZ"].append(hit.GetZ())
    
    data = pd.DataFrame(data_points)
    return data

def MakeTChain():

    # f = r.TFile.Open(args.outputfile, 'recreate')
    # f.cd()
    ch = r.TChain('cbmsim')
    eos = "root://eosuser.cern.ch/"
    basePath = []
    if args.genie:
        if args.inputfile:
            basePath = ["/eos/experiment/sndlhc/MonteCarlo/Neutrinos/Genie/sndlhc_13TeV_down_volTarget_100fb-1_SNDG18_02a_01_000/2/sndLHC.Genie-TGeant4.root"]
            #basePath = ["/eos/experiment/sndlhc/MonteCarlo/Neutrinos/Genie/sndlhc_13TeV_down_volTarget_100fb-1_SNDG18_02a_01_000/1/sndlhc_+volTarget_0.781e16_SNDG18_02a_01_000.0.ghep.root"]
        else:
            for i in range(1, 3):  # Loop over directories 
                file_path = f"/eos/experiment/sndlhc/MonteCarlo/Neutrinos/Genie/sndlhc_13TeV_down_volTarget_100fb-1_SNDG18_02a_01_000/{i}/sndLHC.Genie-TGeant4.root"
                if os.path.exists(file_path):
                    basePath.append(file_path)

        for base in basePath:
            ch.Add(str(base))    

        snd_geo = SndlhcGeo.GeoInterface("/eos/experiment/sndlhc/convertedData/physics/2022/geofile_sndlhc_TI18_V0_2022.root")
        scifi, mufilter = snd_geo.modules['Scifi'], snd_geo.modules['MuFilter']

        return ch, scifi, mufilter

    else:
        if args.muonDIS:
            for i in range(1, 100):  # Loop over directories 
                file_path = f"/eos/experiment/sndlhc/users/dancc/MuonDIS/ecut1.0_z-7_2.5m_Ioni_latelateFLUKA/muonDis_201/{i}/sndLHC.muonDIS-TGeant4-muonDis_201.root"
                if os.path.exists(file_path):
                    basePath.append(file_path)
            
        if args.muonPassing:
            for i in range(1, 27):
                file_path = f"/eos/experiment/sndlhc/MonteCarlo/MuonBackground/muons_up/scoring_2.5/7135377/{i}/sndLHC.Ntuple-TGeant4.root"
                 # file_path = "/eos/experiment/sndlhc/MonteCarlo/MuonBackground/muons_up/scoring_2.5/7135377/1/sndLHC.Ntuple-TGeant4.root"
                basePath.append(file_path)

        if args.neutralhadrons:
            for j in (1, 10):
                file_path = f"/eos/experiment/sndlhc/MonteCarlo/NeutralHadrons/QGSP_BERT_HP_PEN/neutrons/neu_5_10/Ntuples/{j}/sndLHC.PG_2112-TGeant4_20240126.root"
                basePath.append(file_path)

                for i in range(10, 90):
                    file_path = f"/eos/experiment/sndlhc/MonteCarlo/NeutralHadrons/QGSP_BERT_HP_PEN/neutrons/neu_{i}_{i+10}/Ntuples/{j}/sndLHC.PG_2112-TGeant4_20240126.root"
                    if os.path.exists(file_path):
                        basePath.append(file_path)

                file_path = f"/eos/experiment/sndlhc/MonteCarlo/NeutralHadrons/QGSP_BERT_HP_PEN/neutrons/neu_100_150/Ntuples/{j}/sndLHC.PG_2112-TGeant4.root"
                basePath.append(file_path)

                file_path = f"/eos/experiment/sndlhc/MonteCarlo/NeutralHadrons/QGSP_BERT_HP_PEN/neutrons/neu_150_200/Ntuples/{j}/sndLHC.PG_2112-TGeant4.root"
                basePath.append(file_path)

        for base in basePath:
            ch.Add(str(base))    

        return ch

def SaveData(us_data, scifi_data, quark_data):
    scifi_df = pd.DataFrame(scifi_data)
    us_df = pd.DataFrame(us_data)
    quark_df = pd.DataFrame(quark_data)

    quark_df.to_csv(f"{eospath}data_quark_3.csv")
    print(f'Quark csv written to: {eospath}data_quark_3.csv')

    scifi_df.to_csv(f"{eospath}data_scifi_3.csv")
    print(f'Scifi csv written to: {eospath}data_scifi_3.csv')

    us_df.to_csv(f"{eospath}data_us_3.csv")
    print(f'HCAL csv written to: {eospath}data_us_3.csv')

#get numu data from file
if args.genie:
    r.FairTask
    ch, scifi, mufilter = MakeTChain()
    output_file = "/afs/cern.ch/user/t/tismith/sndsw/macro/nu_genie/output"
    if args.inputfile:
        input_file = args.inputfile
        run_convert_script(input_file, output_file)
    else:
        run_convert_script(ch, output_file)
    converted_data = read_convert_script_output(output_file)
    us_data, scifi_data = EventLoop()
    SaveData(us_data, scifi_data, converted_data)

#get data from background file
if args.muonPassing or args.muonDIS or args.neutralhadrons:
    ch = MakeTChain()
    data = get_background_data(ch)
    filename = "PassingMuons_new.csv2"
    # filename = "MuonDIS_100.csv"
    # filename = "NeutralHadrons_largerdataset.csv"

    data.to_csv(f"{eospath}{filename}")
    print(f'csv written to: {eospath}{filename}')

    # muonfile = r.TFile.Open(args.muonPassing, "READ")
    #file = ROOT.TFile.Open("root://eospublic.cern.ch//eos/experiment/sndlhc/MonteCarlo/FLUKA/muons_up/version1/unit30_Nm.root")

    #nt;2 = current cycle, nt;1 = backup cycle
    # ntuple = muonfile.Get("nt;2")
    # ntuple = muonfile.Get("nt;1")

    # df_PassingMuons = getMuonData(ntuple)

    # df_PassingMuons.to_csv(f"{eospath}PassingMuons_otherfile.csv")
    # print(f'Passing Muons csv written to: {eospath}PassingMuons_otherfile.csv')

    # DISfile.Close()

if args.HTCondor:
    ch, scifi, mufilter = MakeTChain()
    us_signal, scifi_signal = GetData(ch, N)
    SaveData(us_data, scifi_data)


