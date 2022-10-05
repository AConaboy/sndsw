import ROOT, os, csv, subprocess, cProfile, io, pstats
from argparse import ArgumentParser
import rootUtils as ut
import AnalysisFunctions as muAna
import SndlhcGeo, SndlhcTracking
from pathlib import Path

deciles=[i/10 for i in range(11)]
v1, v2 = ROOT.TVector3(), ROOT.TVector3()

class MuonSelectionCriteria(ROOT.FairTask):

    def Init(self, options, monitor):
       
        self.M=monitor
        self.options=options
        if self.options.path.find('TI18')>0: self.path='TI18'
        elif self.options.find('H8')>0: self.path='H8'
        run=ROOT.FairRunAna.Instance()
        self.trackTask=run.GetTask('simpleTracking')

        ioman=ROOT.FairRootManager.Instance()
        self.OT=ioman.GetSink().GetOutTree()

        lsOfGlobals=ROOT.gROOT.GetListOfGlobals()
        self.MuFilter=lsOfGlobals.FindObject('MuFilter')
        self.Scifi=lsOfGlobals.FindObject('Scifi')
        self.runNr = str(options.runNumber).zfill(6)
        self.afswork=options.afswork
        self.afsuser=options.afsuser
        self.outpath=self.afswork
        self.EventNumber=-1

        self.systemAndPlanes = {1:2,2:5,3:4}
        #self.systemAndBars={1:7,2:10,3:60}
        #self.systemAndChannels={1:[8,0],2:[6,2],3:[1,0]}
        self.sdict={0:'Scifi',1:'Veto',2:'US',3:'DS'}
        self.zPos=self.M.zPos
        #self.cutdists=muAna.GetCutDistributions(self.runNr, 'yresidual')
        #iterationdict={'zeroth':0, 'ToF':1, 'TW':2}
        #self.iteration=iterationdict[options.mode]

        self.freq=160.316E6
        self.TDC2ns=1E9/self.freq

        #self.largeSiPMmap={0:0 ,1:1 ,3:2 ,4:3 ,6:4 ,7:5}
        #self.verticalBarDict={0:1, 1:3, 2:5, 3:6}
        #verticalPlanes=list(self.verticalBarDict.values())
        #self.xref=42. # To be set as USbarlength/2.
        #self.correctionfunction=lambda ps, qdc, limit: 1/sum( [ ps[i]*qdc**i for i in range(len(ps)) ] ) 
        self.hists=self.M.h

        self.histtypes=['slopes', 'dy', 'nSiPMs']

        slopetitle = 'Track slopes;slope x [rad];slope y [rad]'
        self.hists['slopes'] = ROOT.TH2F('slopes',slopetitle, 300, -1.5, 1.5, 300, -1.5, 1.5)
        for subsystem in (1,2):
            for plane in range(self.systemAndPlanes[subsystem]):
                
                yrestitle=f'y-residual_{10*subsystem+plane};track y - bar y-midpoint [cm];track y-residual;Counts / 1 cm'
                self.hists[f'dy_{10*subsystem+plane}']=ROOT.TH1F(f'dy_{10*subsystem+plane}', yrestitle, 60, -30, 30)

                nSiPMstitle=f'nSiPMs_{10*subsystem+plane};nSiPMs left;nSiPMs right'
                self.hists[f'nSiPMs_{10*subsystem+plane}']=ROOT.TH2F(f'nSiPMs_{10*subsystem+plane}', nSiPMstitle,9,0,9,9,0,9)

                # for bar in range(subsystemAndPlanes[subsystem]):        

    def GetEntries(self):
        return self.eventTree.GetEntries()

    def ExecuteEvent(self, event):
        
        hists=self.hists
        
        ### 1 fired scintillator per muon subsystem (veto+US for TI18, US for H8)
        if self.path=='TI18':tmp=(1,2)
        elif self.path=='H8': tmp=(2,)
        if not muAna.OneHitPerSystem(event.Digi_MuFilterHits, tmp):return 0
        DST0cc=muAna.GetDSH_average(event.Digi_MuFilterHits)
        DST0=self.DST0cut(DST0cc)
        if DST0==0: return 0

        Reco_MuonTracks=self.M.Reco_MuonTracks
        theTrack=Reco_MuonTracks[0]
        fstate=theTrack.getFittedState()
        pos=fstate.getPos()
        mom=fstate.getMom()

        ### Fill slope hist
        slopeX, slopeY = mom.x()/mom.z(), mom.y()/mom.z()
        hists['slopes'].Fill(slopeX, slopeY)
        
        for hit in event.Digi_MuFilterHits:
            detID=hit.GetDetectorID()
            s,p,b=muAna.parseDetID(detID)
            if s==3: continue
            # if not self.yresidual(detID,pos,mom): continue
            
            zEx=self.zPos['MuFilter'][s*10+p]
            lam=(zEx-pos.z())/mom.z()
            Ex=ROOT.TVector3(pos.x()+lam*mom.x(), pos.y()+lam*mom.y(), pos.z()+lam*mom.z())

            xpred=self.xpred(v1,Ex)
            xpred_R=-xpred

            ### Fill y-residual histogram
            self.Fillyresidual(detID,pos,mom,Ex)

            ### Fill nSiPMs histogram
            self.FillnSiPMs(detID, hit)

    def Fillyresidual(self,detID,pos,mom,Ex):
        s,p,b=muAna.parseDetID(detID)
        key=10*s+p
        z=self.zPos['MuFilter'][key]
        lam=(z-pos.z())/mom.z()

        self.MuFilter.GetPosition(detID,v1,v2)
        dy=Ex.y()-v1.y() # Needs to be checked

        self.hists[f'dy_{key}'].Fill(dy)

    def FillnSiPMs(self, detID, hit):
        nLeft, nRight=muAna.GetnFiredSiPMs(hit)
        s,p,b=muAna.parseDetID(detID)
        key=10*s+p
        self.hists[f'nSiPMs_{key}'].Fill(nLeft, nRight)

    def WriteOutHistograms(self):
        
        for h in self.M.h:

                hist=self.M.h[h]
                outpath=f'{self.afswork}/splitfiles/run{self.runNr}/SelectionCriteria/'
                path_obj=Path(outpath)
                path_obj.mkdir(parents=True, exist_ok=True)
                
                outfile=f'SelectionCriteria_{self.options.nStart}.root'
                f=ROOT.TFile.Open(outpath+outfile, 'update')
                f.WriteObject(hist, hist.GetName(), 'kOverwrite')
                f.Close()
        print(f'{len(self.M.h)} histograms saved to {outpath}{outfile}')

    @staticmethod
    def dycut(hist, nsig=1):    
        dymin=hist.GetMean()-nsig*hist.GetStdDev()
        dymax=hist.GetMean()+nsig*hist.GetStdDev()
        return dymin, dymax

    @staticmethod
    def slopecut(mom, slopecut=0.1):
        slopeX, slopeY=mom.x()/mom.z(), mom.y()/mom.z()
        if abs(slopeX)>slopecut or abs(slopeY)>slopecut: return 0
        else: return 1

    def DST0cut(self,DST0):
        if DST0==-999. or DST0==-998. or DST0==-6237.5: return 0
        DST0ns=DST0*6.25
        if DST0ns<13.75 or DST0ns>15.45: return 0
        return DST0ns 
    
    def xpred(self,v1,Ex):
        return ROOT.TMath.Sqrt((v1.x()-Ex.x())**2+(v1.y()-Ex.y())**2+(v1.z()-Ex.z())**2)


