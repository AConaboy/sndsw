import numpy as np
import ROOT
import matplotlib.pyplot as plt
import glob, os
import pandas as pd

base_path = '/eos/experiment/sndlhc/users/aconsnd/simulation/'

fileloc_dict = {
    "neutrino":f"{base_path}neutrino/data/showerprofiles/sndlhc_13TeV_down_volTarget_100fb-1_SNDG18_02a_01_000/",
    "muonDIS":f"{base_path}muonDIS/data/showerprofiles/",
    "passingmuon":f"{base_path}passingmuon/data/showerprofiles/",
    "neutralhadron":f"{base_path}neutralhadron/data/showerprofiles/"
}

filelist_dict = {
    simMode:glob.glob(os.path.join(fileloc_dict[simMode], "*.csv")) for simMode in fileloc_dict.keys()
}

df = pd.concat(
    [pd.read_csv(f).assign(simMode=simMode) for simMode, filelist in filelist_dict.items() for f in filelist], 
    ignore_index=True
)