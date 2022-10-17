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
        self.nStations=options.nStations
        self.outpath=self.afswork
        self.EventNumber=-1

        self.systemAndPlanes = {1:2,2:5,3:4}
        self.sdict={0:'Scifi',1:'Veto',2:'US',3:'DS'}
        self.zPos=self.M.zPos

        self.freq=160.316E6
        self.TDC2ns=1E9/self.freq

        self.hists=self.M.h

        self.histtypes=['slopes', 'dy', 'nSiPMs','nStations']

        self.hists['nStations'] = ROOT.TH1F('nStations', 'Number of stations used in track fit;Number of stations used;Counts', 10, 0, 10)
        self.hists['nStations'].Fill(options.nStations)

        slopetitle = 'Track slopes;slope x [rad];slope y [rad]'
        self.hists[f'slopes_{self.nStations}stations'] = ROOT.TH2F(f'slopes_{self.nStations}stations',slopetitle, 300, -1.5, 1.5, 300, -1.5, 1.5)
        for subsystem in (1,2):
            for plane in range(self.systemAndPlanes[subsystem]):
                
                yrestitle=f'y-residual_{10*subsystem+plane};track y - bar y-midpoint [cm];track y-residual;Counts / 1 cm'
                self.hists[f'dy_{10*subsystem+plane}_{self.nStations}stations']=ROOT.TH1F(f'dy_{10*subsystem+plane}_{self.nStations}stations', yrestitle, 60, -30, 30)

                nSiPMstitle=f'nSiPMs_{10*subsystem+plane};nSiPMs left;nSiPMs right'
                self.hists[f'nSiPMs_{10*subsystem+plane}_{self.nStations}stations']=ROOT.TH2F(f'nSiPMs_{10*subsystem+plane}_{self.nStations}stations', nSiPMstitle,9,0,9,9,0,9)

    def GetEntries(self):
        return self.eventTree.GetEntries()

    def ExecuteEvent(self, event):
        
        hists=self.hists

        if 'slopes' in self.histtypes:
            Reco_MuonTracks=self.M.Reco_MuonTracks
            theTrack=Reco_MuonTracks[0] # For new software update I have to check the unique_ID of the track
            fstate=theTrack.getFittedState()
            pos=fstate.getPos()
            mom=fstate.getMom()

            ### Fill slope hist
            slopeX, slopeY = mom.x()/mom.z(), mom.y()/mom.z()
            hists[f'slopes_{self.nStations}stations'].Fill(slopeX, slopeY)

        # Returns bool for if the track object is in the acceptance for Veto and US respectively.
        InVeto, InUS=(muAna.InAcceptance(pos, mom, i, self.MuFilter, self.zPos) for i in (1,2))
        if not InUS: return 0

        # Use the correct subsystems for 1 bar / plane cut
        if InVeto: tmp=(1,2)
        else: tmp=(2,)
        if not muAna.OneHitPerSystem(event.Digi_MuFilterHits, tmp): return 0

        for hit in event.Digi_MuFilterHits:
            detID=hit.GetDetectorID()
            s,p,b=muAna.parseDetID(detID)
            if s==3: continue

            if 'slopes' in self.histtypes:
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

        self.MuFilter.GetPosition(detID,v1,v2)
        dy=Ex.y()-v1.y() # Needs to be checked
        s,p,b=muAna.parseDetID(detID)
        key=10*s+p
        z=self.zPos['MuFilter'][key]
        lam=(z-pos.z())/mom.z()
        pq = v1-pos
        uCrossv= (v2-v1).Cross(mom)
        doca = pq.Dot(uCrossv)/uCrossv.Mag()

        self.hists[f'dy_{key}_{self.nStations}stations'].Fill(doca)

    def FillnSiPMs(self, detID, hit):
        nLeft, nRight=muAna.GetnFiredSiPMs(hit)
        s,p,b=muAna.parseDetID(detID)
        key=10*s+p
        self.hists[f'nSiPMs_{key}_{self.nStations}stations'].Fill(nLeft, nRight)

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
    
    def xpred(self,v1,Ex):
        return ROOT.TMath.Sqrt((v1.x()-Ex.x())**2+(v1.y()-Ex.y())**2+(v1.z()-Ex.z())**2)


