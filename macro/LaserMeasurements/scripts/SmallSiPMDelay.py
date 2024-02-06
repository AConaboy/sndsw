#!/usr/bin/env python3
import ROOT
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import numpy as np

# Set the default font size for titles
plt.rcParams['figure.titlesize'] = 18
# Set the default font size for titles
plt.rcParams['axes.titlesize'] = 16
# Set the default font size for labels
plt.rcParams['axes.labelsize'] = 16

class SmallSiPMDelays(object):
    def __init__(self, M):
        self.M=M
        self.options=M.options
        if not M.SiPMsize=='small':
            print(f'Selected SiPM is a large SiPM')
            return 
        
        self.delayfoldername='next'
        
    def SetSiPM(SiPM):
        self.M.SetSiPM(SiPM)        
        
    def GetDelays(self):
        self.data={}
        
        for i in range(1,4):
            mode=str(i)
            if not mode in self.data: self.data[mode]={}
            for intensity in self.M.runs[mode]:
                runNumber = self.M.runs[mode][intensity]
                hist=self.GetDelayHist(runNumber)
                if hist==-999:
                    print(f'No hist for intensity {intensity}')
                if not hist or hist==-999: continue
                self.data[mode][intensity] = hist 

    def GetDelayHist(self, runNumber):
        filename = f'{self.M.afswork}plots/run_{str(runNumber).zfill(6)}/plots.root'
        f = ROOT.TFile.Open(filename,'read')
        if not hasattr(f, self.delayfoldername):
            print(f'No folder {self.delayfoldername}')
            return 
        d=f.Get(self.delayfoldername)
        d.cd()
        # histname = f'next_smallSiPM_timestamp{self.M.SiPM}'
        histname = f'next_smallSiPM_timestamp_{self.M.SiPM}'
        if not hasattr(d, histname):
            # print(f'No {histname} hist')
            return -999
        hist = d.Get(histname).Clone() 
        hist.SetDirectory(ROOT.gROOT)
        f.Close()
        return hist
    
    def MakeArrays(self):
        self.data_arrays={}
        for mode in self.data:
            for idx, intensity in enumerate(self.data[mode]):
                if not intensity in self.data_arrays: self.data_arrays[intensity] = np.zeros((30, 50))
                
                h1 = self.data[mode][intensity]
                content = np.array([h1.GetBinContent(xbin) for xbin in range(h1.GetNbinsX())])
                nExpectedSiPMs = int(mode)*8
                self.data_arrays[intensity][nExpectedSiPMs, :] = content 
        
    def PlotSaveDelays(self):
        
        rows=max([len(self.data[i]) for i in self.data])
        fig,axes = plt.subplots(nrows=rows, ncols=1, figsize=(9, rows*9) )

        for idx, intensity in enumerate(self.data_arrays):
            im = axes[idx].imshow(self.data_arrays[intensity], cmap='viridis', origin='lower', aspect='auto')
            cbar = fig.colorbar(im, ax=axes[idx], orientation='vertical', pad=0.1)    
            
            axes[idx].set_title(f'Illuminated SiPMs v small SiPM delay in clock cycles\n{intensity}% on SiPM {self.M.SiPM}')
            axes[idx].set_xlabel('$\Delta$ clock cycles (small - large)')
            axes[idx].set_ylabel('Number of illuminated SiPMs')
            
            axes[idx].set_xticks([i*5 for i in range(50//5)])
            axes[idx].set_yticks([i*8 for i in range(1,4)])
            
            correlationfactor = self.GetCorrelationFactor(intensity)
        
        
        path2fig=f'{self.M.afswork}analysis-plots/SmallSiPMdelays/smallSiPMdelays_SiPM{self.M.SiPM}.png'
        fig.savefig(path2fig)
        print(f'Fig save to {path2fig}')

    def GetCorrelationFactor(self, intensity):
        
        arr = self.data_arrays[intensity]
        dCC_proj, nSiPMs_proj = np.sum(arr, axis=0), np.sum(arr, axis=1)
        
        E_dCC = np.average(np.arange(len(dCC_proj)), weights=dCC_proj)
        E_nSiPMs = np.average(np.arange(len(nSiPMs_proj)), weights=nSiPMs_proj)

        covariance = np.mean(dCC_proj - E_dCC) * np.mean(nSiPMs_proj - E_nSiPMs)
        
        # var_dCC = 

    def PlotDelays(self):
        self.GetDelays()
        self.MakeArrays()
        self.PlotSaveDelays()