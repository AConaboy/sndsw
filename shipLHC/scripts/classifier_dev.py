import numpy as np
import pandas as pd 
import ast, os, csv, shutil

datapath = '/eos/experiment/sndlhc/users/aconsnd/simulation/neutrino/data/nue-extendedreconstruction/'
# cols=['filekey','EventNumber','hasMuon','DSmult0x','DSmult0y','DSmult1x','DSmult1y','USmult3','USmult4','dx3','dy3','dx4','dy4','lambdax3','lambdax4','lambday3','lambday4']
column_names=['filekey', 'EventNumber', 'flav', 'hasMuon', 'fired_planes','interactionWall',
'scifi_median_x','scifi_median_y',
'scifi_residual_x','scifi_residual_y',
'dx0','dx1','dx2','dx3','dx4',
'dy0','dy1','dy2','dy3','dy4',
'ds0','ds1','ds2','ds3','ds4', 'ds_scifi',
'x0','x1','x2','x3','x4',
'y0','y1','y2','y3','y4',
'lambdax0','lambdax1', 'lambdax2','lambdax3', 'lambdax4', 
'lambday0','lambday1', 'lambday2','lambday3', 'lambday4',
'HCAL5barcode', 'm_x','c_x','m_y','c_y'
]

def GetFiles():
    all_files = [i for i in os.listdir(datapath) if all([i.endswith('.csv'), i.find('total')==-1, i.find('FNs')==i.find('FPs')==i.find('TPs')])]
    return all_files

def MakeDf(all_files):
    dfs=[]
    for file in all_files:
        with open(datapath+file, 'r') as f:
            reader=csv.reader(f)
            next(reader)

            # Read lines and apply ast.literal_eval to each line
            # data = [ast.literal_eval(line) for line in reader]
            data = [line for line in reader]
        # print(datapath+file)
        # Step 3: Convert the list of tuples into a pandas DataFrame
        df = pd.DataFrame(data, columns=column_names)
        
        # Step 4: Append the DataFrame to the list
        dfs.append(df)

    combined_df = pd.concat(dfs, ignore_index=True)
    return combined_df

def DfManipulations(df):

    # Some columns have bools to show they couldn't be calculated. I want to change to col type to float and the bools to Nans
    cols_with_bools = ['dx0','dy0','dx1','dy1','dx2','dy2','dx3','dy3','dx4','dy4','lambdax3','lambday3','lambday3','lambday4']

    for col in df.columns: 
        df[col].fillna(0, inplace=True)

    # Calc. the resultant residual 
    for i in range(5):
        print(f'ds for plane {i}')
        df[f'ds{i}'] = np.sqrt(df[f'dx{i}']**2 + df[f'dy{i}']**2)

def SaveDf(df):
    swan_path = '/eos/home-a/aconsnd/SWAN_projects/Data analysis/total_extendedreconstruction.csv'
    df.to_csv(swan_path)
    print(f'Full csv file written to {swan_path}')

def DoAll():
    af=GetFiles()
    df=MakeDf(af)
    # DfManipulations(df)
    SaveDf(df)

def MakeFalseNegPosDirectories():

    # Define the source folders
    source_folders = ['/afs/cern.ch/user/a/aconsnd/Pictures/EventDisplays/FalseNegatives/', '/afs/cern.ch/user/a/aconsnd/Pictures/EventDisplays/FalsePositives/', '/afs/cern.ch/user/a/aconsnd/Pictures/EventDisplays/TruePositives/']
    # Define the destination folder
    destination_folder = '/afs/cern.ch/user/a/aconsnd/Pictures/EventDisplays/'

    # Loop through each source folder
    for source_folder in source_folders:
        # List all files in the current source folder
        for filename in os.listdir(source_folder):
            # Check if the file is a PNG file
            if filename.endswith('.png'):
                # Construct full file path
                source_file = os.path.join(source_folder, filename)
                destination_file = os.path.join(destination_folder, filename)
                
                # Move the file to the destination folder
                shutil.move(source_file, destination_file)
    print(f'Moved all FNs and FPs up 1 directory')

    FNs=[]
    FPs=[]
    TPs=[]
    with open('/eos/experiment/sndlhc/users/aconsnd/simulation/neutrino/data/nue-extendedreconstruction/extendedreconstruction_FNs.csv') as f:
        reader=csv.reader(f)
        next(reader)
        [FNs.append(r) for r in reader]
    with open('/eos/experiment/sndlhc/users/aconsnd/simulation/neutrino/data/nue-extendedreconstruction/extendedreconstruction_FPs.csv') as f:
        reader=csv.reader(f)
        next(reader)
        [FPs.append(r) for r in reader]
    with open('/eos/experiment/sndlhc/users/aconsnd/simulation/neutrino/data/nue-extendedreconstruction/extendedreconstruction_TPs.csv') as f:
        reader=csv.reader(f)
        next(reader)
        [TPs.append(r) for r in reader]    

    eventdisplays = os.listdir('/afs/cern.ch/user/a/aconsnd/Pictures/EventDisplays')
    nuefilter_events = [i for i in eventdisplays if all([i.endswith('png'), i.find('nueFilter')==0])]

    for ev in nuefilter_events:

        filekey, eventNumber = ev.removesuffix('.png').removeprefix('nueFilter_').split('_')
        for FN_key,FN_event in FNs:
            if FN_key==filekey and FN_event==eventNumber:
                os.rename(f'/afs/cern.ch/user/a/aconsnd/Pictures/EventDisplays/nueFilter_{FN_key}_{FN_event}.png', 
                        f'/afs/cern.ch/user/a/aconsnd/Pictures/EventDisplays/FalseNegatives/nueFilter_{FN_key}_{FN_event}.png')
                print(f'Moved event display file {FN_key}, event {FN_event} into FN folder')

        for FP_key,FP_event in FPs:
            if FP_key==filekey and FP_event==eventNumber: 
                os.rename(f'/afs/cern.ch/user/a/aconsnd/Pictures/EventDisplays/nueFilter_{FP_key}_{FP_event}.png', 
                        f'/afs/cern.ch/user/a/aconsnd/Pictures/EventDisplays/FalsePositives/nueFilter_{FP_key}_{FP_event}.png')             
                print(f'Moved event display file {FP_key}, event {FP_event} into FP folder')

        for TP_key,TP_event in TPs:
            if TP_key==filekey and TP_event==eventNumber: 
                os.rename(f'/afs/cern.ch/user/a/aconsnd/Pictures/EventDisplays/nueFilter_{TP_key}_{TP_event}.png',
                        f'/afs/cern.ch/user/a/aconsnd/Pictures/EventDisplays/FalsePositives/nueFilter_{TP_key}_{TP_event}.png')             
                print(f'Moved event display file {TP_key}, event {TP_event} into FP folder')



