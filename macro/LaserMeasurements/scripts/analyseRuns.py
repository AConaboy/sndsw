#!/usr/bin/env python3

import argparse, os, ROOT, csv, json, time
from reverseMapping import reverseChannelMapping
import pandas as pd 
pd.set_option('mode.chained_assignment', None)
import numpy as np
import json
from QDCcalibration import QDCcalibration
from TimeWalk import TimeWalk
from PCBdelays import PCBdelays
from SmallSiPMDelay import SmallSiPMDelays
from AnalysisFunctions import Analysis as LaserAna

parser = argparse.ArgumentParser()
parser.add_argument('--SiPM', dest='SiPM',type=int, default=20)
parser.add_argument('--bar', dest='bar',type=int, default=-1)
parser.add_argument('--PCB', dest='PCB',type=str, default='US')
parser.add_argument('--trigger', dest='trigger',type=str, default='ext')
parser.add_argument('-m','--mode', dest='mode',type=str, default='1')
parser.add_argument('--triggerphase', dest='triggerphase',type=int, default=100)
parser.add_argument('--LaserMeasurements', action='store_true', default=True)
parser.add_argument('--laser-mode', dest='laser_mode', type=str, default='adc')
parser.add_argument('--fitmode', dest='fitmode', type=str, default='linear')
parser.add_argument('--eth', dest='eth', action='store_true')

# Plot making flags
parser.add_argument('--timewalk', action='store_true')
parser.add_argument('--qdccalib', action='store_true')
parser.add_argument('--all-qdccalib', action='store_true')
parser.add_argument('--offset-determination', action='store_true')
parser.add_argument('--make-hists', action='store_true')
parser.add_argument('--pcbdelays', action='store_true')
parser.add_argument('--qdc-intensity-overlay', action='store_true')
parser.add_argument('--timewalk-overlay', action='store_true')
parser.add_argument('--timewalk-overlay-all', action='store_true')
parser.add_argument('--TI18-timewalk-comparison', action='store_true')
parser.add_argument('--TI18-timewalk-comparison-all', action='store_true')
parser.add_argument('--qdc-compare-runtypes', action='store_true')
parser.add_argument('--tw-compare-runtypes', action='store_true')
parser.add_argument('--plot-smallSiPMdelays', action='store_true')

parser.add_argument('--determine-offsets', action='store_true')

options = parser.parse_args()

class LaserRunsData(object):
    def __init__(self,options):
        """ 
        Produce plots for run passed on CL
        """
        self.sndswpath = 'sndsw/macro/LaserMeasurements/'
        
        self.mode=options.mode
        self.path = os.getcwd()
        self.config = self.ReadConfig()
        self.eospath=self.config["eospath"]
        self.datapath=f'{self.eospath}data'
        
        self.SiPM = options.SiPM
        self.TDC2ns = 1E9/160.316E6
        self.options=options
        self.badrunkeywords=self.Badrunkeywords()
        self.GetRunDB() # Makes self.df and self.bars_df

        self.mapping = reverseChannelMapping()
        self.mapping.Init(f'{self.datapath}/run_000000/')
        self.PCB = options.PCB
        if self.PCB=='US': self.TofpetMap = self.mapping.TofpetMap[1]
        elif self.PCB=='DS': self.TofpetMap = self.mapping.TofpetMap[2]
        self.intensityprofile=self.GetIntensityProfile()
        self.h={}
        self.GetSiPMsize()
        self.triggerphase=options.triggerphase

        # Read in calibrated intensities
        self.GetCalibratedIntensities()

        # AnalysisFunctions class for some useful methods
        self.LaserAna = LaserAna(options)

        """
        Dictionary of all runs with keys corresponding 
        to what the laser illuminates: keys= ['SiPM', 1, 2, 3] (bars)
        """
        self.GetRuns()
    
        self.titledict = {'SiPM':'SiPM', '1':'1-bar', '2':'2-bar', '3':'3-bar'}

    def ReadConfig(self):
        # configfile = f"{self.path}/configuration/config.json"
        configfile = f"{self.sndswpath}configuration/config.json"
        with open(configfile, 'r') as f:
            d = json.load(f) 
        return d

    def GetCalibratedIntensities(self):
        calibrated_intensities_file=f"{self.sndswpath}configuration/calibrated_intensities.csv"
        self.calibrated_intensities = pd.read_csv(calibrated_intensities_file, header=0, index_col=False)     

    def SetSiPM(self, SiPM):
        self.SiPM = SiPM
        x=self.GetSiPMsize()
        if x==-999: return x
        self.GetRuns()
        
    def SetMode(self, mode):
        self.mode=mode

    def GetSiPMsize(self):
        condition=self.df['SiPM number']==self.SiPM 
        fvals=list(self.df.loc[condition, 'SiPM size'])
        if len(fvals)==0:
            print(f'No valid runs using SiPM {self.SiPM}')
            return -999
        self.SiPMsize = list(self.df.loc[condition, 'SiPM size'])[0]
    
    def GetSiPMRuns(self):
        # Runs with the desired SiPM
        condition1 = self.df['SiPM number']==self.SiPM 
        condition2 = self.df['Trigger phase (cc)']==self.triggerphase
        
        self.SiPMruns={}
        # Loop through intensity values, find runs at those intensities for these SiPMs and append
        for intensity in self.intensityprofile[f'{self.PCB}-{self.SiPMsize}']:
            condition3 = self.df['Laser intensity (%)']==intensity
            fval = self.df.loc[condition1 & condition2 & condition3, 'RunNumber']
            if len(fval)==0: 
                print(f'SiPM {self.SiPM} missing intensity {intensity}') 
                continue
            self.SiPMruns[intensity]=list(fval)[-1]
            
    def string_to_list(self, s):
        return [int(x) for x in s.split(',')]
        
    def generate_exp_sipms(self, bars):
        return [list(range(int(1 + (b - 1) * 8), b * 8 + 1)) for b in bars]
        
    def GetRuns(self):
        """
        Find SiPM, 1-bar, 2-bar and 3-bar runs for all intensities
        """

        self.bar, self.SiPM_in_group = self.LaserAna.SiPM2BarAndPosition(self.SiPM)
        bar_df = self.bars_df[self.bars_df['bar'].str.find(str(self.bar)) != -1]
        # contains_bar = lambda x : any([int(i) == self.bar for i in x.split(',')])
        contains_bar = lambda x : self.bar in x 

        bar_df_copy = bar_df.copy()
        bar_df.loc[:,"has_bar"] = bar_df_copy["bar"].apply(contains_bar)

        self.bar_df = bar_df[bar_df['has_bar']]

        self.runs={}
        self.GetSiPMRuns()
        self.runs['SiPM'] = self.SiPMruns

        # Bar runs are taken with large SiPM intensity profile
        for intensity in self.intensityprofile[f'{self.PCB}-large']: 
            condition1 = bar_df['Laser intensity (%)']==intensity
            for i in (1,2,3):

                if not str(i) in self.runs: self.runs[str(i)]={}

                # condition2 = self.bar_df['bar'].str.split(',').apply(len) == i
                condition2 = self.bar_df['bar'].apply(len) == i
                fvals = self.bar_df.loc[ condition1 & condition2, "RunNumber"]

                if len(fvals)==0: continue
                self.runs[str(i)][intensity]=list(fvals)[-1]

    def CheckForIntensities(self):
        if not hasattr(self, runs): self.GetRuns()
        measured_intensities=list(self.SiPMruns.keys())
        for intensity in self.intensityprofile:
            if not intensity in measured_intensities: print(f'{intensity} data not found')

    def GetIntensityProfile(self):
        with open(f'{self.sndswpath}configuration/intensityprofile.csv', 'r') as f:
            d={i[0]:i[1:] for i in csv.reader(f)}
        for i in d: d[i] = [float(j) for j in d[i]]
        return d
    
    def GetUSSiPMpositions(self):
        # with open(f'{path}configuration/US-SiPMpositions.txt', 'r') as f:
        df=pd.read_csv(f'{self.path}/configuration/US-SiPMpositions.txt', sep=' ', names=['element','name', 'library', 'package', 'value', 'x', 'y', 'locked', 'smashed', 'rot'])
        df.drop(columns=['element', 'library', 'package', 'value', 'locked', 'smashed', 'rot'])
        self.USSiPMpositions=df
            
    def GetHistogramValue(self, runNr, mode):
        
        filename = f'{self.eospath}plots/run_{str(runNr).zfill(6)}/plots.root'
        
        if not os.path.exists(filename): 
            print(f'filename does not exist: {filename}\nMaking plots...')
            self.AnalyseRun(runNr)
         
        isSmall = self.LaserAna.IsSmallSiPMchannel(self.SiPM_in_group-1) 
   
        f=ROOT.TFile.Open(filename, 'read')
        if not hasattr(f, 'qdcvtimestamp'):
            print(f'{mode} plots not produced for run {runNr}')
            f.Close()
            return
        tdir=f.Get('qdcvtimestamp')
        tdir.cd()
        if not hasattr(tdir, f'qdcvtimestamp_SiPM{self.SiPM}'):
            f.Close()
            return   
        
        histogram = tdir.Get(f'qdcvtimestamp_SiPM{self.SiPM}')
        if options.eth:
            if histogram.GetEntries() < 500: 
                print(f'SiPM {self.SiPM}, run {runNr} fails entry threshold of 500')
                f.Close()
                return            
                
        xproj=histogram.ProjectionX() # QDC distribution
        modevalue = xproj.GetBinCenter(xproj.GetMaximumBin()) # Find most popular QDC value
        stddev=xproj.GetStdDev()
        binlow, binhigh = xproj.FindBin(modevalue-stddev), xproj.FindBin(modevalue+stddev)
        if mode=='QDC': 
            trunc_entries = sum([xproj.GetBinContent(i) for i in range(binlow, binhigh)])
            if trunc_entries==0: return
            stderror = modevalue / ROOT.TMath.Sqrt(trunc_entries) *1/0.67 # approx qdc dist as Gaussian
            return modevalue, stderror
        elif mode=='timestamp':
            y_sig = histogram.ProjectionY('y_sig', binlow, binhigh)
            timestamp_modevalue = y_sig.GetBinCenter(y_sig.GetMaximumBin()) # Find most popular QDC value
            return timestamp_modevalue, y_sig.GetMeanError()

    def AnalyseRun(self, runNr):
        # options_subset
        options.runNumber=runNr
        options.make_hists=True
        from analyseRun import LaserData
        ld=LaserData(options)
        ld.EventLoop()
        ld.WriteOutHistograms()
        print("Made hists using LaserData")        

    def GetLastRun(self):
        runfiles=os.listdir(self.datapath+'/data/')
        lastrun=max([int(i.split('_')[-1]) for i in runfiles if len(i.split('.'))==1])
        return str(lastrun).zfill(6)

    def Badrunkeywords(self):
        with open(f'{self.sndswpath}configuration/badrunkeywords.csv', 'r') as f:
            d=list(csv.reader(f))[0]
        return d         

    def GetRunDB(self):
        configfiles=[f for f in os.listdir(f'{self.sndswpath}configuration/') if f.split('.')[-1]=='csv']
        for f in configfiles: 
            if f.find(options.PCB)>0: filename=f

        df = pd.read_csv(f'{self.sndswpath}configuration/{filename}', sep=',', header=0)

        # Remove rows that start with a comment
        df=df.dropna(subset=['SiPM number', "Laser intensity (%)", 'Optics set up'])

        # Require matching of the k:v pairs in self.path/configuration/config.json
        for k in self.config:
            if not k in df.columns: continue
            df=df.drop(df[df[k] != self.config[k]].index)
        
        # Remove runs by keyword flag in Run status column
        for badword in self.badrunkeywords:
            df=df.drop(df[df['Run status'].str.lower() == badword].index)
            
        # Write run number as type integer
        df['RunNumber'] = df['RunNumber'].astype(int)
        
        # Drop runs with multi bars just for the SiPM runDB
        self.df=df.drop(df[df['SiPM number'].str.find(',')!=-1].index)
        
        # Write SiPM number as type integer
        self.df['SiPM number'] = self.df['SiPM number'].astype(int)

        # Copy df as bars_df when there is atleast 1 group illuminated
        self.bars_df=df.drop(df[df['SiPM size']!='bar'].index)
        
        # Drop bar runs from SiPM DB
        self.df=self.df.drop(self.df[self.df['SiPM size'] == 'bar' ].index)
        
        # Rename column for multi-SiPM DB
        self.bars_df = self.bars_df.rename(columns={'SiPM number': 'bar'})
        
        # Write bars column as type list
        self.bars_df['bar'] = self.bars_df['bar'].apply(self.string_to_list)
        
        # Create column of expected SiPMs to fire for bars_df
        self.bars_df['exp. SiPMs'] = self.bars_df['bar'].apply(lambda x: self.generate_exp_sipms(x))

        
    def WriteOutGoodRunNrs(self):
        unique_runNrs = pd.concat([self.df['RunNumber'], self.bars_df['RunNumber']]).unique()
        
        outfile = f'{self.sndswpath}configuration/GoodRuns.txt'
        with open(outfile, 'w') as f:
            for r in unique_runNrs:
                f.write(f'{r}\n')
                
    def WriteOutBarRunNrs(self):
        
        selected_columns = self.bars_df[['RunNumber', 'exp. SiPMs']]
        
        outfile = f'{self.sndswpath}configuration/SingleBarRuns.txt'
        
        selected_columns.to_csv(outfile, header=None, index=False)
        
    def FindSiPMsize(self, SiPM):
        return self.df.loc[self.df['SiPM number']==SiPM, 'SiPM size'].values[0]

def timewalk(modes=[options.mode]):
    for mode in modes: 
        tw.SetMode(mode)
        tw.TimeWalk()  

def QDCcalib():
    qdccalib.MPVprogression()    

def all_qdccalib():
    for SiPM in range(1, 81):
        print(f"Plotting MPV progression for SiPM {SiPM}")
        qdccalib.SetSiPM(SiPM)
        qdccalib.MPVprogression()  
    qdccalib.PlotOffsets()      

def offset_determination():
    for SiPM in range(1,81):
        print(f"Determining QDC offset for SiPM {SiPM}")
        qdccalib.SetSiPM(SiPM)
        qdccalib.DetermineQDCoffset()
    qdccalib.PlotOffsets()
    qdccalib.WriteOutOffsets()
    
    # Write out offsets
    offsets_filename = f'{lds.afswork}Results/Offsets/qdcoffsets_{lds.titledict[lds.mode]}.csv'
    offsets_df = pd.DataFrame(qdccalib.linear_params, columns=['SiPM number', 'Linear fit'])
    offsets_df.to_csv(offsets_filename, index=False)    

def qdc_intensity_overlay():
    if options.bar==-1: bar=[i for i in range(1,11)]
    else: bar=[options.bar]
    for b in bar: 
        print(f'Overlaying QDC-intensity plots for bar {b}')
        qdccalib.OverlaySiPMs(bar=b)    

def qdc_intensity_overlay():
    if options.bar==-1: bar=[i for i in range(1,11)]
    else: bar=[options.bar]
    for b in bar: 
        print(f'Overlaying QDC-intensity plots for bar {b}')
        qdccalib.OverlaySiPMs(bar=b)    

def qdc_compare_runtypes():
    if options.bar==-1: bar=[2,3,4,5,6,7]
    else: bar=[options.bar]
    for b in bar: 
        SiPMs_in_bar = [i for i in range((b-1)*8+1, b*8+1)]
        for SiPM in SiPMs_in_bar: 
            lds.SetSiPM(SiPM)
            print(f'Overlaying QDC-intensity plots for SiPM {SiPM}')
            qdccalib.CompareRunTypes()
def tw_compare_runtypes():
    if options.bar==-1: bar=[2,3,4,5,6,7]
    else: bar=[options.bar]
    for b in bar: 
        SiPMs_in_bar = [i for i in range((b-1)*8+1, b*8+1)]
        for SiPM in SiPMs_in_bar: 
            lds.SetSiPM(SiPM)
            print(f'Overlaying QDC-intensity plots for SiPM {SiPM}')
            tw.CompareRunTypes()            

def timewalk_overlay(modes=[options.mode]):
    
    bar=lds.LaserAna.SiPM2BarAndPosition(lds.SiPM)[0]
    for mode in modes:
        tw.SetMode(mode)
        tw.OverlaySiPMs(bar)

def timewalk_overlay_all():

    for i in ('SiPM', '1', '2', '3'):
        if not i =='3': bars=[i for i in range(1,11)]
        else: bars=[2,3,4,5,6,7]
        for b in bars:
            tw.SetMode(i)
            tw.OverlaySiPMs(b)

def TI18_timewalk_comparison():
    tw.TI18_laser_comparison()
    
def TI18_timewalk_comparison_allmodes():
    for i in ('SiPM', '1', '2', '3'):
        tw.SetMode(i)
        TI18_timewalk_comparison_all()

def TI18_timewalk_comparison_all():
    if lds.mode!='3':SiPMrange=range(1,81)
    else: SiPMrange=range(9,49)
    for SiPM in SiPMrange:
        print(f'SiPM {SiPM}')
        tw.SetSiPM(SiPM)
        if lds.SiPMsize=='small': continue
        tw.TI18_laser_comparison()
    tw.PlotSystemTI18LaserDifference()
        

def pcbdelays():
    delays.GetData()
    delays.PlotDelays()  
    
def plot_smallSiPMdelays():
    smallSiPMdelays.PlotDelays()

tic=time.perf_counter()
lds = LaserRunsData(options)
tw, qdccalib, delays, smallSiPMdelays = TimeWalk(lds), QDCcalibration(lds), PCBdelays(lds), SmallSiPMDelays(lds)
toc=time.perf_counter()
print(f'Made all classes in: {round(toc-tic,2)}s')

if options.timewalk: 
    timewalk()

if options.qdccalib:
    QDCcalib()

elif options.all_qdccalib:
    all_qdccalib()

if options.offset_determination:
    offset_determination()

if options.qdc_intensity_overlay:
    qdc_intensity_overlay()

if options.qdc_compare_runtypes:
    qdc_compare_runtypes()
if options.tw_compare_runtypes:
    tw_compare_runtypes()    

if options.timewalk_overlay:
    timewalk_overlay()
if options.timewalk_overlay_all:
    timewalk_overlay_all()
        
if options.TI18_timewalk_comparison:
    TI18_timewalk_comparison()
    
if options.TI18_timewalk_comparison_all:
    TI18_timewalk_comparison_all() 


if options.pcbdelays:
    pcbdelays()

if options.plot_smallSiPMdelays:
    plot_smallSiPMdelays()