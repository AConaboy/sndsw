#!/usr/bin/env python3

import argparse,os,ROOT,uproot,time
from datetime import datetime
from reverseMapping import reverseChannelMapping
import pandas as pd
pd.set_option('mode.chained_assignment', None)
import awkward as ak
import numpy as np
from ast import literal_eval

path='/eos/experiment/sndlhc/raw_data/commissioning/US_tests_LabLausanne_2023/'
afsoutpath='/afs/cern.ch/work/a/aconsnd/LaserMeasurements/plots/'

parser = argparse.ArgumentParser()
parser.add_argument('--runNumber', '-r', dest='runNumber',type=int, default=-1)
parser.add_argument('--PCB', dest='PCB',type=str, default='US')
parser.add_argument('--make-hists', action='store_true')
parser.add_argument('--make-dt-hists', action='store_true')
options = parser.parse_args()

ROOT.gStyle.SetOptStat(111111)

class LaserData(object):
    def __init__(self,options):
        """ 
        Produce plots for run passed on CL
        """

        self.options=options
        if options.runNumber>0: self.runNr = str(options.runNumber).zfill(6)
        else: 
            self.runNr = GetLastRun()
        self.datapath=f'{path}data/run_{self.runNr}/'

        """ 
        Load in dataframe of expected SiPMs for each run
        If run in self.exp_SiPMs_df then get the expected small and large SiPMs
        also set options.make-dt-hists to True
        """
        self.exp_SiPMs_df = pd.read_csv(f'/afs/cern.ch/user/a/aconsnd/LaserMeasurements/configuration/BarRuns.txt', names=['RunNumber', 'Exp. SiPMs'])
        if options.runNumber in self.exp_SiPMs_df['RunNumber'].values: 
            print(f'Bar run. Making dt hists also.')
            self.GetExpectedSiPMs()
            options.make_dt_hists=True

        self.mapping = reverseChannelMapping()
        self.mapping.Init(self.datapath)
        if options.PCB=='US': self.TofpetMap = self.mapping.TofpetMap[1]
        if options.PCB=='DS': self.TofpetMap = self.mapping.TofpetMap[2]
        self.h={}
        
        self.GetData()

    def SiPM2BarAndPosition(self, SiPM):
        # Pass the SiPM number on the PCB 
        # returns the bar number and SiPM number within the bar
        barNumber = (SiPM-1)//8+1
        SiPMNumber = (SiPM-1)%8+1
        if SiPMNumber==0: SiPMNumber==8
        return barNumber, SiPMNumber

    def IsSmallSiPM(self, SiPMNumber):
        if SiPMNumber in (3,6): return True
        else: return False
        
    def string_to_list(self, i):
        return [int(x) for x in i if x.isdigit()]     
        
    def GetExpectedSiPMs(self):
    
        self.exp_SiPMs_df['Exp. SiPMs'] = self.exp_SiPMs_df['Exp. SiPMs'].apply(literal_eval)
        self.exp_SiPMgroups = self.exp_SiPMs_df.loc[self.exp_SiPMs_df['RunNumber']==int(self.runNr), 'Exp. SiPMs'].array[0]
        self.exp_SiPMs = [item for si in self.exp_SiPMgroups for item in si]
        self.exp_smallSiPMs = [ [i[2], i[5]] for i in self.exp_SiPMgroups]
        self.exp_smallSiPMs = [item for si in self.exp_smallSiPMs for item in si] 
        self.exp_largeSiPMs = [i for i in self.exp_SiPMs if i not in self.exp_smallSiPMs]       

    def GetData(self, cut=False):

        tic=time.perf_counter()
        
        filename=f"{self.datapath}data_0000.root"
        self.df = pd.DataFrame(uproot.open(filename)['event_data'].arrays(library='np'))    

        """
        Find timestamp peak to use as a QCD filter.
        Just load in the timestamp hist for the SiPM, take a rough range of mode +/- 2 std dev 
        and use this as a cut on QDC to remove noise.
        """
        if cut:
            mode_timestamp = self.Find_timestamp_peak()

            # # Select a window around the peak_timestamp +/- 200 ps
            # mask = (self.df['timestamp'].apply(lambda x: all(value >= mode_timestamp-0.2 and value <= mode_timestamp+0.2 for value in x)))

            # self.timing_condition_df = self.df[mask]
            # self.timing_condition_df['SiPM number'] = self.timing_condition_df.apply(
            #     lambda row: self.GetSiPMNumber(row['tofpet_id'], row['tofpet_channel']),axis=1)
        else:
            self.timing_condition_df = self.df
            self.timing_condition_df['SiPM number'] = self.timing_condition_df.apply(
                lambda row: self.GetSiPMNumber(row), axis=1)
                        
        toc=time.perf_counter()
        print(f'Time elapsed: {toc-tic:0.4f} seconds')
        
    def GetSiPMNumber(self, row):
        return np.array([self.TofpetMap[t_id*1000 + t_ch] for t_id,t_ch in zip(row['tofpet_id'], row['tofpet_channel'])])
        
    def EventLoop(self):

        """
        Histograms to make:
        1. qdc for all fired SiPMs
        2. timestamps for all fired SiPMs 
        3. 2D plot of qdc for all SiPMs
        3. Correlation of timestamp and saturation
        """    
        
        self.FillHists()
        if options.make_dt_hists:   self.dt_hists()

    def FillHists(self):

        for index, row in self.timing_condition_df.iterrows():
            
            data = zip(row['SiPM number'], row['value'], row['timestamp'], row['value_saturation'], row['t_coarse'], row['v_coarse'])

            PCBqdc_histname=f'QDC_PCB'
            if not PCBqdc_histname in self.h:
                title = 'QDC distribution for all SiPMs on PCB;QDC [a.u];SiPM number'
                self.h[PCBqdc_histname] = ROOT.TH2F(PCBqdc_histname, title, 290, -20, 270, 80, 1, 81)
            # [self.h[PCBqdc_histname].Fill(qdc, SiPM_number) for SiPM_number, qdc in data]
    
            for SiPM, qdc, timestamp, value_saturation, t_coarse, v_coarse in data: 
                bar,sipm_in_group = self.SiPM2BarAndPosition(SiPM)
                isSmall = self.IsSmallSiPM(sipm_in_group)
                
                self.h[PCBqdc_histname].Fill(qdc, SiPM)
                
                # if isSmall:
                qdcvtimestamp_histname=f'qdcvtimestamp_SiPM{SiPM}'
                if not qdcvtimestamp_histname in self.h:
                    title = 'QDC v timestamp for SiPM '+str(SiPM)+';QDC [a.u];Timestamp [cc]'
                    self.h[qdcvtimestamp_histname] = ROOT.TH2F(qdcvtimestamp_histname, title, 290, -20, 270, 500, 0, 5)
                self.h[qdcvtimestamp_histname].Fill(qdc,timestamp)                
                
                saturation_histname=f'saturation_SiPM{SiPM}'
                if not saturation_histname in self.h:
                    title = 'Saturation distribution for SiPM '+str(SiPM)+';QDC saturation (%);Counts'
                    self.h[saturation_histname] = ROOT.TH1F(saturation_histname, title, 100, 0, 1)
                self.h[saturation_histname].Fill(value_saturation)
                
                vcoarse_histname=f'vcoarse_SiPM{SiPM}'
                if not vcoarse_histname in self.h:
                    title = 'v_coarse distribution for SiPM '+str(SiPM)+';v_coarse;Counts'
                    self.h[vcoarse_histname] = ROOT.TH1F(vcoarse_histname, title, 60, 0, 60)
                self.h[vcoarse_histname].Fill(v_coarse)
                
                tcoarse_histname=f'tcoarse_SiPM{SiPM}'
                if not tcoarse_histname in self.h:
                    title = 't_coarse distribution for SiPM '+str(SiPM)+';t_coarse;Counts'
                    self.h[tcoarse_histname] = ROOT.TH1F(tcoarse_histname, title, 60, 0, 60)
                self.h[tcoarse_histname].Fill(t_coarse) 
                
    def dt_hists(self):
        
        # Used for comparing the time that small and large SiPMs fire
        for self.index, row in self.timing_condition_df.iterrows():
            
            self.SetEventTimeStamp()
            # data = zip(row['SiPM number'], row['value'], row['timestamp'], row['value_saturation'])
            
            fractionExpSiPMs_histname=f'frac_SiPMs'
            if not fractionExpSiPMs_histname in self.h:
                title = 'Fraction of expected SiPMs that fire;N_{fired} / N_{expected};Counts'
                self.h[fractionExpSiPMs_histname] = ROOT.TH1F(fractionExpSiPMs_histname, title, len(self.exp_SiPMs), 0, 1)
            
            n_large, n_small = self.GetFractionExpectedSiPMs(row['SiPM number'])
            frac_n = (n_large+n_small) / len(self.exp_SiPMs)
            self.h[fractionExpSiPMs_histname].Fill(frac_n) 
            
            """
            Only makes sense to look at when the small SiPMs are firing relative to the current event, if the current event contains 
            a certain fraction of the expected large SiPMs! 
            
            To find this I first look at the distribution of the fraction of the expected fired SiPMs
            """
            
            next_smallSiPM_tcoarse_histname=f'next_smallSiPM_tcoarse'
            if not next_smallSiPM_tcoarse_histname in self.h:
                title = 't_coarse of delayed small SiPM hit;t_coarse;Counts'
                self.h[next_smallSiPM_tcoarse_histname] = ROOT.TH1I(next_smallSiPM_tcoarse_histname, title, 5, 0, 5)
                
            next_smallSiPM_vcoarse_histname=f'next_smallSiPM_vcoarse'
            if not next_smallSiPM_vcoarse_histname in self.h:
                title = 'v_coarse of delayed small SiPM hit;v_coarse;Counts'
                self.h[next_smallSiPM_vcoarse_histname] = ROOT.TH1F(next_smallSiPM_vcoarse_histname, title, 50, 0, 50) 
            
            if n_large >= 0.5*len(self.exp_SiPMs):
                nClockCycles = self.GetNextSmallSiPMHit(row['SiPM number']) # SiPM: [nCC, QDC]
                for i in nClockCycles.items():
                    SiPM, data = i
                    clockCycle, smallQDC, v_coarse, t_coarse = data 
                    self.h[next_smallSiPM_tcoarse_histname].Fill(t_coarse)
                    self.h[next_smallSiPM_vcoarse_histname].Fill(v_coarse)
                    
                    next_smallSiPM_timestamp_histname=f'next_smallSiPM_timestamp_{SiPM}'
                    if not next_smallSiPM_timestamp_histname in self.h:
                        title = f'Clock cycles until expected small SiPM {SiPM} next fires;Clock cycles [n];Counts'
                        self.h[next_smallSiPM_timestamp_histname] = ROOT.TH1I(next_smallSiPM_timestamp_histname, title, 50, 0, 50)
                    self.h[next_smallSiPM_timestamp_histname].Fill(clockCycle)

                    # What is happening for the peak at 5 clock cycles later?
                    if clockCycle != 0: 
                        LateSmallSiPM_histname=f'LateSmallSiPM_{SiPM}'
                        if not LateSmallSiPM_histname in self.h:
                            title = 'QDC of delayed small SiPM hits;QDC [a.u];Counts'
                            self.h[LateSmallSiPM_histname] = ROOT.TH1F(LateSmallSiPM_histname, title, 50, -20, 30)
                        
                        self.h[LateSmallSiPM_histname].Fill(smallQDC)
                
                # Fill overflow bin for each exp. small SiPM not found in 100 events
                if len(nClockCycles) < len(self.exp_smallSiPMs):
                    for i in range( len(self.exp_smallSiPMs) - len(nClockCycles) ): self.h[next_smallSiPM_timestamp_histname].Fill(self.h[next_smallSiPM_timestamp_histname].GetNbinsX()+1)

    def SetEventTimeStamp(self):
        self.evt_timestamp = self.timing_condition_df.iloc[self.index]['evt_timestamp']
            
    def GetFractionExpectedSiPMs(self, SiPMs):
        n_exp = len(self.exp_SiPMs) 
        n_large=0
        n_small=0
        for i in SiPMs: 
            if i in self.exp_largeSiPMs: n_large+=1
            if i in self.exp_smallSiPMs: n_small+=1
        return n_large, n_small
    
    def GetNextSmallSiPMHit(self, SiPMs):
        
        """
        Steps to find the event timestamp of the next events when the expected small SiPMs fire
        1a. Check if the small SiPMs fire in the same event as when >= 0.5* exp. large SiPMs fire
        1b. If so, fill the delta t hist with 0 for each of the small SiPM that fire in the same event
        1c. If all the expected small SiPMs fire in the same event, return.
        
        2a. Loop over the subsequent events of the run (changing to next 100 events) and look for the remaining small SiPMs
        2b. If a small SiPM is found and has not already been accounted for in the original event or a previous event, get that event timestamp
        2c. If all expected small SiPMs are accounted for, return
        
        3. Fill delta t histogram with the difference in event time between the large SiPM event, and the event where the expected small SiPMs fire
        
        """
        
        clockCycles = {} 
        
        # Check whether all the expected small SiPMs fire in the event: if so, return 
        for idx, SiPM in enumerate(SiPMs): 
            if SiPM in self.exp_smallSiPMs and SiPM not in clockCycles: 
                smallSiPMQDC = self.timing_condition_df['value'].iloc[self.index][idx]
                smallSiPMvcoarse = self.timing_condition_df['v_coarse'].iloc[self.index][idx]
                smallSiPMtcoarse = self.timing_condition_df['t_coarse'].iloc[self.index][idx]
                clockCycles[SiPM]=[0, smallSiPMQDC, smallSiPMvcoarse, smallSiPMtcoarse]

        if len(clockCycles)==len(self.exp_smallSiPMs): return clockCycles

        else:
            index=self.index
            
            remaining_events = self.timing_condition_df.shape[0] - index - 1
            for evt_number in range( min(100, remaining_events) ):

                row = self.timing_condition_df.iloc[self.index+evt_number]
                for idx, SiPM in enumerate(row['SiPM number']): 
                    if SiPM in self.exp_smallSiPMs and not SiPM in clockCycles: 
                        
                        delta_evt_timestamp = row['evt_timestamp'] - self.evt_timestamp 
                        smallSiPMQDC = row['value'][idx]
                        smallSiPMvcoarse = row['v_coarse'][idx]
                        smallSiPMtcoarse = row['t_coarse'][idx]
                    
                        clockCycles[SiPM] = [delta_evt_timestamp, smallSiPMQDC, smallSiPMvcoarse, smallSiPMtcoarse]

                if len(clockCycles) == len(self.exp_smallSiPMs): break 
        
        return clockCycles
            
    def Find_timestamp_peak(self):
        exploded_df = self.df['timestamp'].explode()
        mode_timestamp = exploded_df.mode()[0]
        return mode_timestamp

    def WriteOutHistograms(self):

        plotspath=f'{afsoutpath}run_{self.runNr}/'
        if not os.path.exists(plotspath):os.mkdir(plotspath)

        timinghistkeys = ['frac_SiPMs', 'next_smallSiPM_timestamp']

        outfilename=f'{plotspath}plots.root'
        outfile=ROOT.TFile.Open(outfilename, 'recreate')
        for hname in self.h:
            if hname not in timinghistkeys:
                key=hname.split('_')[0]
                if not hasattr(outfile, key): outfile.mkdir(key)

            else: 
                key='timing'
                if not hasattr(outfile, key): outfile.mkdir(key)

            tdir=outfile.Get(key)
            tdir.cd()

            hist=self.h[hname]
            hist.Write(hname, 2)

        outfile.Close()
        print(f'Histograms written to {outfilename}')

def GetLastRun():
    runfiles=os.listdir(f'{path}data/')
    lastrun=max([int(i.split('_')[-1]) for i in runfiles if len(i.split('.'))==1])
    return str(lastrun).zfill(6)
    
if options.make_hists or options.make_dt_hists:
    ld=LaserData(options)
    ld.EventLoop()
    ld.WriteOutHistograms()