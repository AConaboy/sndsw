import csv, os
import pandas as pd 
from matplotlib import pyplot as plt 
import xgboost as xgb 
import numpy as np

from sklearn.metrics import make_scorer, accuracy_score, f1_score, recall_score, precision_score, classification_report
from sklearn.model_selection import train_test_split, GridSearchCV, learning_curve, StratifiedKFold
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay, roc_curve, auc, precision_recall_curve, matthews_corrcoef
from sklearn.inspection import permutation_importance

base_path = '/eos/experiment/sndlhc/users/aconsnd/simulation/neutrino/data/showerprofiles/sndlhc_13TeV_down_volTarget_100fb-1_SNDG18_02a_01_000/'

# weight is a constant for the neutrino simulation
columns = ['hasMuon','interactionWall','lambdax0','lambdax1','lambdax2','lambdax3','lambdax4','lambday0','lambday1','lambday2','lambday3','lambday4','relQDC0','relQDC1','relQDC2','relQDC3','relQDC4','planeQDC0','planeQDC1','planeQDC2','planeQDC3','planeQDC4','hasTrack','vetoCut','avgSFchannel','avgDSchannel','SF1','SF2','consec2SFplanes','allHCALplanes']

def GetFiles():
    csv_files = [f for f in os.listdir(base_path) if f.endswith('.csv')]
    return csv_files 

def GetData(files):

    dataframes = []

    # Decile-based progress tracking
    decile_markers = {int(np.ceil(len(files) * p)) for p in np.linspace(0.1, 1.0, 10)}

    # Loop through each CSV file
    for idx,csv_file in enumerate(files):
        print(f'{idx}/{len(files)}')
        # Read the CSV file
        file_path = os.path.join(base_path, csv_file)
        df = pd.read_csv(file_path)[columns]
        
        #'vetoCut','avgSFchannel','avgDSchannel','SF1','SF2','consec2SFplanes','allHCALplanes'
        df_filtered = df[(df['vetoCut']) & (df['allHCALplanes']) & (df['consec2SFplanes']) & (df['SF1']) & (df['SF2'])]
        
        # Append the filtered DataFrame to the list
        dataframes.append(df_filtered)

        # Print heartbeat at decile checkpoints
        if idx in decile_markers:
            print(f"Progress: {int((idx / len(files)) * 100)}% ({idx}/{len(files)} files processed)")        

    # Concatenate all filtered DataFrames into a single DataFrame
    combined_df = pd.concat(dataframes, ignore_index=True)    

    return combined_df

### Optimising model ### 
def BestParameters(mode='optimise'):
    
    cv_strategy = StratifiedKFold(n_splits=10, shuffle=True, random_state=42)
    # Make custom scorer ensuring that pos_label=0
    f1_scorer = make_scorer(custom_f1_scorer)

    if mode=='optimise':

        # Define the parameter grid for GridSearchCV
        param_grid = {
            'n_estimators': [100, 200],
            'learning_rate': [0.01, 0.05],
            'max_depth': [2, 3],
            'min_child_weight': [1, 5],  # Helps to control overfitting
            'scale_pos_weight': [scale_pos_weight]            
        }
        
        # Instance model
        bdt_model = xgb.XGBClassifier(eval_metric='logloss', random_state=42)
        
        # Instantiate GridSearchCV
        grid_search = GridSearchCV(
            estimator=bdt_model, 
            param_grid=param_grid, 
            cv=cv_strategy, 
            scoring=f1_scorer, 
            verbose=1)

        # Fit the model using grid search
        grid_search.fit(X, y)

        # Get the best parameters
        best_params = grid_search.best_params_
        features = list(X.columns)
    elif mode=='load':
        with open('/afs/cern.ch/user/a/aconsnd/best_numuBDTparams.json') as json_file:
            best_params = json.load(json_file)
        features = best_params.pop('features')
    print(f"Best Parameters: {best_params}")
    return best_params, features

def custom_f1_scorer(y_true, y_pred):
    return f1_score(y_true, y_pred, pos_label=0) 
def custom_recall_scorer(y_true, y_pred):
    return recall_score(y_true, y_pred, pos_label=0, zero_division=0)
def custom_precision_scorer(y_true, y_pred):
    return precision_score(y_true, y_pred, pos_label=0, zero_division=0) 
def custom_accuracy_scorer(y_true, y_pred):
    return accuracy_score(y_true, y_pred, pos_label=0)  

def GetXY(df):
    X = df.drop(['hasMuon', 'relQDC0', 'relQDC1', 'relQDC2', 'relQDC3', 'relQDC4', 'planeQDC0', 'planeQDC1', 'planeQDC2', 'planeQDC3', 'planeQDC4', 'vetoCut', 'avgSFchannel', 'avgDSchannel', 'SF1', 'SF2',       'consec2SFplanes', 'allHCALplanes'], axis=1)    
    y=df['hasMuon']
    X.columns
    return X, y   