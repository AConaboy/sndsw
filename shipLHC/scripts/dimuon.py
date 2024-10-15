#!/usr/bin/env python
import ROOT, os

class Dimuon(object):

    def __init__(self, options, tw):
       
        # self.state='corrected'
        # options.state='corrected'
        self.options=options
        self.tw=tw
        self.runNr = tw.runNr 
        self.muAna = tw.muAna
        self.timealignment=tw.timealignment

        ### If no time-walk correction run is provided. Set the default correction run depending on time alignment of the data set

        self.afswork=tw.afswork
        self.outpath=tw.outpath

        self.subsystemdict={1:'Veto', 2:'US', 3:'DS'}
        self.nchs={1:224, 2:800}

        self.systemAndPlanes = {1:2,2:5,3:7}
        self.systemAndBars={1:7,2:10,3:60}
        self.systemAndChannels={1:[8,0],2:[6,2],3:[1,0]}
        self.sdict={0:'Scifi',1:'Veto',2:'US',3:'DS'}
        self.zPos=tw.zPos
        self.cutdists=self.muAna.GetCutDistributions(self.runNr, ('dy', 'timingdiscriminant'))

        self.freq=160.316E6
        self.TDC2ns=1E9/self.freq

        self.xrefs=tw.xrefs
        self.sides=('left', 'right')

        self.hists=tw.hists 
        
    def FillHists(self, hits):
        """
        Events passed here have 2 DS tracks. The tds0 is not determined yet and no requirement
        on the red.chi2 of either track is placed yet. 
        """
        
        # Get DS event t0. All fired DS bars are used to determine this.
        x=self.muAna.GetDSHaverage(hits) # Now returns in ns
        if x==-999.: 
            print(f'Event {self.tw.M.EventNumber} has no DS horizontal hits with 2 fired SiPMs')
            return
        self.TDS0, firedDSHbars=x
        if not 'TDS0' in self.hists:
           self.hists['TDS0']=ROOT.TH1F('TDS0','Average time of DS horizontal bars;DS horizontal average time [ns];Counts', 200, 0, 50)
        self.hists['TDS0'].Fill(self.TDS0) 
        
        if not 'nFiredDSbars' in self.hists:
            self.hists['nFiredDSbars'] = ROOT.TH1F('nFiredDSbars','Fired DS bars;Number of fired bars in the DS;Counts', 40, 0, 40)
        self.hists['nFiredDSbars'].Fill(firedDSHbars)
        
        
        
        
        