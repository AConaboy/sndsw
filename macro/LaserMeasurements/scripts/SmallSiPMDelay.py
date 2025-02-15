#!/usr/bin/env python3
import ROOT, os
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import numpy as np

# Set the default font size for titles
plt.rcParams['figure.titlesize'] = 18
# Set the default font size for titles
plt.rcParams['axes.titlesize'] = 16
# Set the default font size for labels
plt.rcParams['axes.labelsize'] = 16
fontsize=16

class SmallSiPMDelays(object):
    def __init__(self, M):
        self.M=M
        self.options=M.options
        
        self.delayfoldername = 'next'
        self.intensity=60 # the intensity point to find delay with
        
    def SetSiPM(self, SiPM):
        self.M.SetSiPM(SiPM)
        
    def GetDelays(self):
        self.data={}
        for SiPM in range(1, 81):
            # only studying small SiPMs
            if self.M.FindSiPMsize(SiPM)=='large': continue
            self.SetSiPM(SiPM)
        
            runNumber = self.M.runs[self.M.mode][self.intensity]
            delay = self.GetTimingDelay(runNumber)
            if not delay or delay==-999: continue
            self.data[SiPM] = delay

    def GetTimingDelay(self, runNumber):
        filename = f'/afs/cern.ch/work/a/aconsnd/LaserMeasurements/plots/run_{str(runNumber).zfill(6)}/plots.root'
        f = ROOT.TFile.Open(filename,'read')
        if not hasattr(f, self.delayfoldername):
            print(f'No folder {self.delayfoldername}')
            return -999
        d=f.Get(self.delayfoldername)
        d.cd()
        # histname = f'next_smallSiPM_timestamp{self.M.SiPM}'
        histname = f'next_smallSiPM_timestamp_{self.M.SiPM}'
        if not hasattr(d, histname):
            print(f'No {histname} hist')
            return -999
        hist = d.Get(histname).Clone() 
        hist.SetDirectory(ROOT.gROOT)
        delay = self.FindDelay(hist)
        return delay

    def FindDelay(self, hist):

        qdc_proj = hist.ProjectionY()
        if qdc_proj.GetMean()==qdc_proj.GetStdDev()==0: return -999
        modal_qdc_bin = qdc_proj.GetMaximumBin()
        modal_qdc = qdc_proj.GetBinCenter(modal_qdc_bin)

        timing_proj_at_qdcmode = hist.ProjectionX('timing_proj_at_qdcmode', modal_qdc_bin, modal_qdc_bin)

        timing_peak = timing_proj_at_qdcmode.GetBinCenter(timing_proj_at_qdcmode.GetMaximumBin())

        return timing_peak
    
    # def MakeArrays(self):
    #     self.data_arrays={}
    #     for mode in self.data:
    #         for idx, intensity in enumerate(self.data[mode]):
    #             if not intensity in self.data_arrays: self.data_arrays[intensity] = np.zeros((30, 50))
                
    #             h1 = self.data[mode][intensity]
    #             content = np.array([h1.GetBinContent(xbin) for xbin in range(h1.GetNbinsX())])
    #             nExpectedSiPMs = int(mode)*8
    #             self.data_arrays[intensity][nExpectedSiPMs, :] = content 

    def PlotDelays(self):
        fig,ax = plt.subplots(figsize=(5,5))
        ax.set_title(f'Delay in clock cycles until expected small SiPM fires\nfor {self.intensity}% laser intensity, {self.M.titledict[self.M.mode]} illumination',fontsize=fontsize)
        ax.tick_params(axis='both', which='major', labelsize=fontsize)

        smallSiPMs = list(self.data.keys())
        # ax.set_xticks(smallSiPMs)
        # ax.set_xticklabels(smallSiPMs, fontsize=10)
        [ax.axvline(SiPM, linestyle='--', color='grey', alpha=0.5) for SiPM in smallSiPMs]

        ax.set_xlabel('SiPM number',fontsize=fontsize )
        ax.set_ylabel('Delay [cc]',fontsize=fontsize)
        # ax.set_ylim(0, 20)
        ax.scatter(self.data.keys(), self.data.values())
        
        plt.tight_layout()
        d=f'{self.M.sndswpath}analysis-plots/smallSiPMdelays/'
        os.makedirs(d, exist_ok=True)
        plotlocation = d+f'{self.M.titledict[self.M.mode]}-{self.intensity}int.png'

        plt.savefig(plotlocation, bbox_inches='tight')
        print(f'Plot saved to {plotlocation}')
    # def PlotSaveDelays(self):
        
        # rows=max([len(self.data[i]) for i in self.data])
        # fig,axes = plt.subplots(nrows=rows, ncols=1, figsize=(9, rows*9) )

        # for idx, intensity in enumerate(self.data_arrays):
        #     im = axes[idx].imshow(self.data_arrays[intensity], cmap='viridis', origin='lower', aspect='auto')
        #     cbar = fig.colorbar(im, ax=axes[idx], orientation='vertical', pad=0.1)    
            
        #     axes[idx].set_title(f'Illuminated SiPMs v small SiPM delay in clock cycles\n{intensity}% on SiPM {self.M.SiPM}')
        #     axes[idx].set_xlabel('$\Delta$ clock cycles (small - large)')
        #     axes[idx].set_ylabel('Number of illuminated SiPMs')
            
        #     axes[idx].set_xticks([i*5 for i in range(50//5)])
        #     axes[idx].set_yticks([i*8 for i in range(1,4)])
            
        #     correlationfactor = self.GetCorrelationFactor(intensity)
        
        
        # path2fig=f'{self.M.afswork}analysis-plots/SmallSiPMdelays/smallSiPMdelays_SiPM{self.M.SiPM}.png'
        # fig.savefig(path2fig)
        # print(f'Fig save to {path2fig}')

    def GetCorrelationFactor(self, intensity):
        
        arr = self.data_arrays[intensity]
        dCC_proj, nSiPMs_proj = np.sum(arr, axis=0), np.sum(arr, axis=1)
        
        E_dCC = np.average(np.arange(len(dCC_proj)), weights=dCC_proj)
        E_nSiPMs = np.average(np.arange(len(nSiPMs_proj)), weights=nSiPMs_proj)

        covariance = np.mean(dCC_proj - E_dCC) * np.mean(nSiPMs_proj - E_nSiPMs)
        