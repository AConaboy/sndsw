import csv, os
import pandas as pd 
from matplotlib import pyplot as plt 
import xgboost as xgb 

base_path = '/eos/experiment/sndlhc/users/aconsnd/simulation/neutrino/data/showerprofiles/sndlhc_13TeV_down_volTarget_100fb-1_SNDG18_02a_01_000/'

# weight is a constant for the neutrino simulation
columns = ['hasMuon','interactionWall','lambdax0','lambdax1','lambdax2','lambdax3','lambdax4','lambday0','lambday1','lambday2','lambday3','lambday4','relQDC0','relQDC1','relQDC2','relQDC3','relQDC4','planeQDC0','planeQDC1','planeQDC2','planeQDC3','planeQDC4','hasTrack','vetoCut','avgSFchannel','avgDSchannel','SF1','SF2','consec2SFplanes','allHCALplanes']

def GetFiles():
    csv_files = [f for f in os.listdir(base_path) if f.endswith('.csv')]
    return csv_files 

def GetData(files):

    dataframes = []
    # Loop through each CSV file
    for idx,csv_file in enumerate(files):
        print(f'idx/{len(files)}')
        # Read the CSV file
        file_path = os.path.join(base_path, csv_file)
        df = pd.read_csv(file_path)[columns]
        
        # Remove rows where 'vetoCut' is "False"
        df_filtered = df[(df['vetoCut']) & (~df['hasMuon']) & (df['allHCALplanes'])]
        
        # Append the filtered DataFrame to the list
        dataframes.append(df_filtered)

    # Concatenate all filtered DataFrames into a single DataFrame
    combined_df = pd.concat(dataframes, ignore_index=True)    

    return combined_df