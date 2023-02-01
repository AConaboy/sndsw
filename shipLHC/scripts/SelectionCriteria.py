import ROOT, os, csv, subprocess, cProfile, io, pstats
from argparse import ArgumentParser
import rootUtils as ut
import AnalysisFunctions as muAna
import SndlhcGeo, SndlhcTracking
from pathlib import Path
from array import array

deciles=[i/10 for i in range(11)]
A, B, locA, locB = ROOT.TVector3(), ROOT.TVector3(), ROOT.TVector3(), ROOT.TVector3()


class MuonSelectionCriteria(ROOT.FairTask):

    def Init(self, options, monitor):
       
        self.M=monitor
        self.options=options
        if self.options.path.find('TI18')>0: self.path='TI18'
        elif self.options.path.find('H8')>0: self.path='H8'
        run=ROOT.FairRunAna.Instance()
        self.trackTask=run.GetTask('simpleTracking')

        ioman=ROOT.FairRootManager.Instance()
        self.OT=ioman.GetSink().GetOutTree()

        lsOfGlobals=ROOT.gROOT.GetListOfGlobals()
        self.MuFilter=lsOfGlobals.FindObject('MuFilter')
        self.Scifi=lsOfGlobals.FindObject('Scifi')
        self.nav=ROOT.gGeoManager.GetCurrentNavigator()
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

        # self.subnamedict={1:'Veto', 2:'US'}
        # for i in (1, 2):
        #     slopetitle = f'{self.subnamedict[i]} track slopes;slope x [rad];slope y [rad]'
        #     self.hists[f'{self.subnamedict[i]}_slopes_{self.nStations}stations'] = ROOT.TH2F(f'{self.subnamedict[i]}_slopes_{self.nStations}stations',slopetitle, 300, -1.5, 1.5, 300, -1.5, 1.5)
        slopetitle = f'Track slopes;slope x [rad];slope y [rad]'
        self.hists[f'slopes_{self.nStations}stations'] = ROOT.TH2F(f'slopes_{self.nStations}stations',slopetitle, 300, -1.5, 1.5, 300, -1.5, 1.5)        

        timingdistitle="DSH average - US1 average;T^{DS}_{0} - #bar{t_{US1}};Counts"
        self.hists['timingdiscriminant']=ROOT.TH1F('timingdiscriminant',timingdistitle, 500,-25,50)

        NPlanestitle='Fraction of expected planes firing with only 1 scintillator;N_{fired}/N_{expected} [dimensionless];Counts'
        self.hists['NPlanes']=ROOT.TH1F('NPlanes', NPlanestitle, 10, 0., 1.)

        self.hists['DST0']=ROOT.TH1F('DST0','DSH average', 100, 0, 25)

        for subsystem in (1,2):
            for plane in range(self.systemAndPlanes[subsystem]):
                
                yrestitle=f'y-residual_{10*subsystem+plane};track y - bar y-midpoint [cm];track y-residual;Counts / 1 cm'
                self.hists[f'dy_{10*subsystem+plane}_{self.nStations}stations']=ROOT.TH1F(f'dy_{10*subsystem+plane}_{self.nStations}stations', yrestitle, 60, -30, 30)

                nSiPMstitle=f'nSiPMs_{10*subsystem+plane};nSiPMs left;nSiPMs right'
                self.hists[f'nSiPMs_{10*subsystem+plane}_{self.nStations}stations']=ROOT.TH2F(f'nSiPMs_{10*subsystem+plane}_{self.nStations}stations', nSiPMstitle,9,0,9,9,0,9)

        title='#chi_{#nu}^{2} of DS track fitting with '+str(options.nStations)+' projections per plane;#chi_{#nu}^{2} [dimensionless];Counts'
        self.hists['trackfitting']=ROOT.TH2F('trackfitting', title, 100,0.0,100.,10,-0.5,9.5)

    def GetEntries(self):
        return self.eventTree.GetEntries()

    def ExecuteEvent(self, event):
        
        hists=self.hists
        hits=event.Digi_MuFilterHits

        tracks={}
        Reco_MuonTracks = self.M.Reco_MuonTracks
        inVeto, inDS=False, False
        for i,track in enumerate(Reco_MuonTracks):
            if track.GetUniqueID()==1: 
                inVeto=True
                tracks[1]=Reco_MuonTracks[i]
            if track.GetUniqueID()==3:
                inDS=True
                tracks[2]=Reco_MuonTracks[i]
        if not inDS: return

        fstate=tracks[2].getFittedState()
        fitStatus=tracks[2].getFitStatus()
        if not fitStatus.isFitConverged(): return
        # fstates={i:tracks[i].getFittedState() for i in tracks}
        pos=fstate.getPos()
        # posvectors={i:fstates[i].getPos() for i in fstates}
        mom=fstate.getMom()
        # momvectors={i:fstates[i].getMom() for i in fstates}

        # Use the correct subsystems for 1 bar / plane cut
        if inVeto: tmp=(1,2,3)
        else: tmp=(2,3)
        # Nfired=muAna.OneHitPerSystem(hits, tmp, Nfired=True)
        # hists['NPlanes'].Fill(Nfired)
        if not muAna.OneHitPerSystem(hits, tmp): return        

        DST0cc=muAna.GetDSHaverage(hits)
        DST0=DST0cc*6.25
        self.hists['DST0'].Fill(DST0)

        if 'slopes' in self.histtypes:
            ### Fill slope hist
            # slopes={i:(momvectors[i].x()/momvectors[i].z(), momvectors[i].y()/momvectors[i].z()) for i in momvectors}
            # [hists[f'{self.subnamedict[i]}slopes_{self.nStations}stations'].Fill(*slopes[i]) for i in slopes] 
            slopes=mom.x()/mom.z(), mom.y()/mom.z()
            hists[f'slopes_{self.nStations}stations'].Fill(*slopes)

        # Returns bool for if the track object is in the acceptance for Veto and US respectively.

        timingdiscriminant=self.GetTimingDiscriminant(hits,DST0)
        hists['timingdiscriminant'].Fill(timingdiscriminant)

        trackchi2NDF=fitStatus.getChi2()/fitStatus.getNdf()+1E-10
        hists['trackfitting'].Fill(trackchi2NDF, fitStatus.getNdf())

        for hit in hits:

            detID=hit.GetDetectorID()
            s,p,b=muAna.parseDetID(detID)
            if s==3: continue

            if 'slopes' in self.histtypes:
                zEx=self.zPos['MuFilter'][s*10+p]
                # lam=(zEx-posvectors[s].z())/momvectors[s].z()
                lam=(zEx-pos.z())/mom.z()
                # Ex=ROOT.TVector3(posvectors[s].x()+lam*momvectors[s].x(), posvectors[s].y()+lam*momvectors[s].y(), posvectors[s].z()+lam*momvectors[s].z())
                Ex=ROOT.TVector3(pos.x()+lam*mom.x(), pos.y()+lam*mom.y(), pos.z()+lam*mom.z())
                
                # Fill scintillator positions into A & B
                self.MuFilter.GetPosition(detID,A,B)

                
                # Find the resultant distance between +ve x-side of scintillator and track 
                xpred=self.xpred(A,Ex)
                xpred_R=-xpred

                ### Fill y-residual histogram
                # self.Fillyresidual(detID,posvectors[s],momvectors[s],Ex)
                self.Fillyresidual(detID,pos,mom,Ex)

            ### Fill nSiPMs histogram
            self.FillnSiPMs(detID, hit)

    def Fillyresidual(self,detID,pos,mom,Ex):

        # dy=Ex.y()-A.y() # Needs to be checked
        # s,p,b=muAna.parseDetID(detID)
        # key=10*s+p
        # z=self.zPos['MuFilter'][key]
        # lam=(z-pos.z())/mom.z()
        # pq = A-pos
        # uCrossv= (B-A).Cross(mom)
        # doca = pq.Dot(uCrossv)/uCrossv.Mag()

        self.MuFilter.GetPosition(detID, A, B)
        self.MuFilter.GetLocalPosition(detID, locA, locB)
        s,p,b=muAna.parseDetID(detID)
        key=10*s+p
        z=self.zPos['MuFilter'][key]
        lam=(z-pos.z())/mom.z()
        
        xEx, yEx, zEx=pos.x()+lam*mom.x(), pos.y()+lam*mom.y(), pos.z()+lam*mom.z()
        Ex=array('d', [xEx, yEx, zEx])
        locEx=array('d', [0,0,0])
        self.nav.MasterToLocal(Ex, locEx)
        locPos=0.5*(locA+locB)
        # Add on half the 
        barY=self.MuFilter.GetConfParF('MuFilter/UpstreamBarY') if s==2 else self.MuFilter.GetConfParF('MuFilter/DownstreamBarY')
        dy=locPos[1]-locEx[1]+barY*0.5 

        self.hists[f'dy_{key}_{self.nStations}stations'].Fill(dy)

    def FillnSiPMs(self, detID, hit):
        nLeft, nRight=muAna.GetnFiredSiPMs(hit)
        s,p,b=muAna.parseDetID(detID)
        key=10*s+p
        self.hists[f'nSiPMs_{key}_{self.nStations}stations'].Fill(nLeft, nRight)

    def GetTimingDiscriminant(self, hits,DST0):
        for hit in hits:
            detID=hit.GetDetectorID()
            s,p,b=muAna.parseDetID(detID)
            if s==2 and p==0:
                US1hit=hit
                break
        averageUS1TDC=muAna.GetAverageTDC(US1hit)
        averageUS1time=averageUS1TDC*6.25
        return DST0-averageUS1time

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
    
    def xpred(self,A,Ex):
        return ROOT.TMath.Sqrt((A.x()-Ex.x())**2+(A.y()-Ex.y())**2+(A.z()-Ex.z())**2)


