#!/usr/bin/env python3

import matplotlib.pyplot as plt

class PCBdelays(object):       
    
    def __init__(self, M):
        self.options = M.options
        self.M = M 
        self.intensity_colour='C1'
        self.qdc_colour='C0'
        self.colours = ['blue', 'green', 'red', 'cyan', 'magenta', 'yellow', 'tab:orange', 'tab:purple']        

        # self.smallcondition1 = self.M.df['SiPM size']=='small'
        # self.smallcondition2 = self.M.df['Laser intensity (%)']==self.M.intensityprofile[f'{self.M.PCB}-small'][-1]
        # self.largecondition1 = self.M.df['SiPM size']=='large'
        # self.largecondition2 = self.M.df['Laser intensity (%)']==self.M.intensityprofile[f'{self.M.PCB}-large'][-1]

        # self.conditions = {"small":[smallcondition1, smallcondition2], "large":[largecondition1, largecondition2]}

        self.maxintensities = {'small':self.M.intensityprofile['US-small'][-1], 'large':self.M.intensityprofile['US-large'][-1]}

    ### Get highest intensity runs for all SiPMs available. 
    ### 

    def GetData(self):
        
        if self.M.PCB=='US': nSiPMs=80
        else: nSiPMs=60

        self.data={}

        for SiPM in range(1,nSiPMs+1):

            # Update LaserRunsData.runs['SiPM']
            print(f'Getting high intensity timestamp for SiPM {SiPM}')
            self.M.SetSiPM(SiPM)

            if self.M.LaserAna.IsSmallSiPMchannel(SiPM-1): runNr = self.M.runs['SiPM'][self.maxintensities['small']]
            else: runNr = self.M.runs['SiPM'][self.maxintensities['large']]

            timestamp = self.M.GetHistogramValue(runNr, 'timestamp')
            
            if not timestamp: 
                print(f'Run {runNr} has no SiPM {self.M.SiPM} plots made')
                continue
            self.data[SiPM] = [i*self.M.TDC2ns for i in timestamp]

    def PlotDelays(self):
        timestamps, timestamps_errors = zip(*self.data.values())
        fig, ax = plt.subplots()

        bars = ax.bar(self.data.keys(), timestamps, yerr=timestamps_errors, capsize=5, color='orange', edgecolor='black')
        # plt.bar(self.data.keys(), self.data.values())

        ax.set_xlabel('SiPM number')
        ax.set_ylabel('Timestamp [ns]')
        ax.set_title('Timestamp recorded for SiPMs with highest laser intensity')
        ax.xaxis.grid(True)

        plt.savefig(f'{self.M.sndswpath}analysis-plots/pcb-delays/pcb-delays.png')

        # Make histogram for PCB delays
        # plt.bar()

