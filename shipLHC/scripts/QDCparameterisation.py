#!/usr/bin/env python
import ROOT, os
import numpy as np
from itertools import combinations

class QDCparameterisation(object):

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
        self.sigmatds0=0.263 # ns 
    
        # self.MuFilter = self.tw.MuFilter
        self.A, self.B = ROOT.TVector3(), ROOT.TVector3()

        self.X, self.Y = 10, 5
        SiPMcombinations_tmp = {'left':list(combinations([0,1,3,4,6,7], 2)), 
                                'right':list(combinations([8,9,11,12,14,15], 2))
                                }
        self.SiPMcombinations = []
        [self.SiPMcombinations.extend(i) for i in SiPMcombinations_tmp.values()]

        self.dict_combinations = {self.SiPMcombinations[idx]:idx for idx in range(len(self.SiPMcombinations))}

        self.shape=(self.X, self.Y, len(self.SiPMcombinations), self.options.nEvents)
        self.MakeDataDict()

        self.USbarX, self.USbarY = self.tw.MuFilter.GetConfParF('MuFilter/UpstreamBarX'), self.tw.MuFilter.GetConfParF('MuFilter/UpstreamBarX')
        self.Xbin, self.Ybin = self.USbarX/self.X, self.USbarY/self.Y

    def MakeDataDict(self):
        self.data = {}
        for p in range(5):
            for b in range(10):
                detID = int(f'2{p}00{b}')
                self.data[detID] = np.zeros(self.shape)

    def Execute(self, hit):
        
        event = self.tw.M.EventNumber - self.options.nStart

        detID = hit.GetDetectorID()
        s,p,b = self.muAna.parseDetID(detID)
        
        zEx=self.tw.zPos['MuFilter'][s*10+p]
        lam=(zEx-self.tw.pos.z())/self.tw.mom.z()
        Ex=ROOT.TVector3(self.tw.pos.x()+lam*self.tw.mom.x(), self.tw.pos.y()+lam*self.tw.mom.y(), self.tw.pos.z()+lam*self.tw.mom.z())

        x_bin, y_bin = self.FindCell(detID, Ex.x(), Ex.y())

        qdcs = hit.GetAllSignals()
        rel_qdcs = self.GetRelativeQDCs(qdcs)

        arr = self.data[detID]

        for comb in rel_qdcs:
            relQDC = rel_qdcs[comb]
            arr[x_bin, y_bin, self.dict_combinations[comb], event]=relQDC
 
    def GetRelativeQDCs(self, qdcs):

        qdcs_dict = {SiPM:qdc for SiPM,qdc in qdcs if not self.muAna.IsSmallSiPMchannel(SiPM)}
        side_wise = {"left":
                            {i:qdcs_dict[i] for i in qdcs_dict if i<8},
                    "right":
                            {i:qdcs_dict[i] for i in qdcs_dict if i>7}
                    }
        comb_dict = {'left':{}, 'right':{}}

        for side in side_wise:
            d=side_wise[side]
            for k1, k2 in combinations(d.keys(), 2):
                relQDC = d[k1]/d[k2]
                comb_dict[side][(k1, k2)] = relQDC

        res={}
        [res.update(sd) for sd in comb_dict.values()]
        return res

        
        


        
        # res = {x:qdcs_dict[x[0]]/qdcs_dict[x[1]] for x in qdcs_dict}

        return res

    def FindCell(self, detID, xEx, yEx):

        self.tw.MuFilter.GetPosition(detID, self.A, self.B)

        y_base = 0.5*(self.A.y()+self.B.y()) - self.USbarY/2
        y_max = 0.5*(self.A.y()+self.B.y()) + self.USbarY/2
        dy = yEx-y_base 
        y_bin = int(dy/self.Ybin)

        # x_base = 0.5*(self.A.x()+self.B.y()) - self.USbarY/2
        x_base = self.A.x()
        x_max = self.B.x()+self.B.x()
        dx = xEx-x_base
        x_bin = int(dx/self.Xbin)

        return x_bin, y_bin

