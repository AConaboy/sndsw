import ROOT, os, csv, subprocess, cProfile, io, pstats
from argparse import ArgumentParser

from AnalysisFunctions import Analysis
from pathlib import Path
from array import array

deciles=[i/10 for i in range(11)]
A, B, locA, locB = ROOT.TVector3(), ROOT.TVector3(), ROOT.TVector3(), ROOT.TVector3()

ROOT.gStyle.SetTitleXOffset(0.7)
ROOT.gStyle.SetTitleYOffset(0.6)
# ROOT.gStyle.SetTitleXSize(0.06)
# ROOT.gStyle.SetTitleYSize(0.06)
# ROOT.gStyle.SetTitleSize(0.12)
for i in ('x', 'y', 'z','t'): ROOT.gStyle.SetTitleSize(0.04, i)

class SystemAlignment(ROOT.FairTask):

    def Init(self, options, monitor):
       
        self.state='corrected'
        options.state='corrected'
        
        self.muAna=Analysis(options)
        self.M=monitor
        self.options=options
        if self.options.path.find('commissioning/TI18')>0:
            self.outpath=options.afswork+'-commissioning/'
            self.path='TI18'
        elif self.options.path.find('physics/2022')>0: 
            self.outpath=options.afswork+'-physics2022/'
            self.path='TI18'
        elif self.options.find('H8')>0: 
            self.outpath=options.afswork+'-H8/'
            self.path='H8'

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
        self.EventNumber=-1
        self.subsystemdict={1:'Veto', 2:'US', 3:'DS'}
        self.nchs={1:224, 2:800}

        self.systemAndPlanes = {1:2,2:5,3:7}
        self.systemAndBars={1:7,2:10,3:60}
        self.systemAndChannels={1:[8,0],2:[6,2],3:[1,0]}
        self.sdict={0:'Scifi',1:'Veto',2:'US',3:'DS'}
        self.zPos=self.M.zPos
        self.cutdists=self.muAna.GetCutDistributions(self.runNr, ('dy', 'timingdiscriminant'))
        self.timealignment=self.muAna.GetTimeAlignmentType(runNr=self.runNr)

        self.freq=160.316E6
        self.TDC2ns=1E9/self.freq

        self.largeSiPMmap={0:0 ,1:1 ,3:2 ,4:3 ,6:4 ,7:5}
        self.verticalBarDict={0:1, 1:3, 2:5, 3:6}
        self.xref=42. # To be set as USbarlength/2.
        
        # self.correctionfunction=lambda ps, qdc: 1/sum( [ ps[i]*qdc**i for i in range(len(ps)) ] ) 
        self.getaverage=lambda d, key, i:sum( [(d[key][i]) for k in range(2) ])/len(d)        

        self.hists=self.M.h
        self.ChannelsWithNoParams=[]
        
        self.CorrectionType=4
        self.correctionparams = lambda ps : [y for x,y in enumerate(ps) if x%2==0]
        self.correctionfunction = lambda ps, qdc : ps[3]*(qdc-ps[0])/( ps[1] + ps[2]*(qdc-ps[0])*(qdc-ps[0]) ) 

        ### If no time-walk correction run is provided. Set the default correction run depending on time alignment of the data set
        if options.TWCorrectionRun==None:
            if self.timealignment=='old': self.TWCorrectionRun=str(5097).zfill(6)
            elif self.timealignment=='new': self.TWCorrectionRun=str(5408).zfill(6)
        ### Check that the time-walk correction run selected has the same time alignment has the data set
        else:
            if self.muAna.GetTimeAlignmentType(runNr=options.TWCorrectionRun) != self.timealignment:
                if self.timealignment=='old': self.TWCorrectionRun=str(5097).zfill(6)
                elif self.timealignment=='new': self.TWCorrectionRun=str(5408).zfill(6)
            else: self.TWCorrectionRun=str(options.TWCorrectionRun).zfill(6)        
        
        # self.cscintvalues=self.muAna.Makecscintdict(self.TWCorrectionRun, state='corrected')
        # self.twparameters=self.muAna.MakeTWCorrectionDict(self.TWCorrectionRun)

    def GetEntries(self):
        return self.eventTree.GetEntries()

    def ExecuteEvent(self, event):

        tracks={1:[], 3:[]}
        Reco_MuonTracks=self.M.Reco_MuonTracks
        inVeto, inDS=False, False
        for i,track in enumerate(Reco_MuonTracks):
            if any([not track.getFitStatus().isFitConverged(), track.getFitStatus().getNdf()==0]): continue
            if track.GetUniqueID()==1: 
                inVeto=True
                tracks[1].append(Reco_MuonTracks[i])
            if track.GetUniqueID()==3: 
                inDS=True
                tracks[3].append(Reco_MuonTracks[i])
        if not inDS: return
        
        if len(tracks[3])==1: dstrack= tracks[3][0]
        else: 
            tmp={i:i.getFitStatus().getChi2()/i.getFitStatus().getNdf() for i in tracks[3]}
            dstrack=tmp[ min(tmp) ]
        
        fitStatus=dstrack.getFitStatus()
        fstate=dstrack.getFittedState()
        
        self.pos=fstate.getPos()
        self.mom=fstate.getMom()
        self.trackchi2NDF=fitStatus.getChi2()/fitStatus.getNdf()+1E-10

        # Get DS event t0
        TDS0=self.muAna.GetDSHaverage(self.MuFilter, event.Digi_MuFilterHits) # Now returns in ns
        if not TDS0: return
        if TDS0==-999.: print(f'Event {self.M.EventNumber} has no DS horizontal hits with 2 fired SiPMs')
        if not 'TDS0' in self.hists:
           self.hists['TDS0']=ROOT.TH1F('TDS0','Average time of DS horizontal bars;DS horizontal average time [ns];Counts', 200, 0, 50)
        self.hists['TDS0'].Fill(TDS0)

        # if self.options.timingdiscriminantcut:
        td = self.GetTimingDiscriminant() # Require that US1 TDC average is less than the DSH TDC average to ensure forward travelling track
        if not td: return
        print(f'{td}')
        if not self.TimingDiscriminantCut(td): return
        
        ### Slope cut
        if not SystemAlignment.slopecut(self.mom): return

        for hit in event.Digi_MuFilterHits:
            nLeft, nRight=self.muAna.GetnFiredSiPMs(hit)
            
            detID=hit.GetDetectorID()
            s,p,b=self.muAna.parseDetID(detID)
            
            # Only investigate veto hits if there is a Scifi track!
            if not inVeto and s==1: continue

            # This function also calls MuFilter.GetPosition with the detID so it's important even for the DS, where the 
            # y-residual cut is not applied.
            self.MuFilter.GetPosition(detID,A,B)
            self.MuFilter.GetLocalPosition(detID, locA, locB)
            if not self.yresidual3(detID): continue

            zEx=self.zPos['MuFilter'][s*10+p]
            lam=(zEx-self.pos.z())/self.mom.z()
            Ex=ROOT.TVector3(self.pos.x()+lam*self.mom.x(), self.pos.y()+lam*self.mom.y(), self.pos.z()+lam*self.mom.z())

            self.GetDistanceToSiPM(A,Ex)

            channels_t=hit.GetAllTimes()
            channels_qdc=hit.GetAllSignals()

            if s!=3: self.FillBarTimeHists(hit)

    def FillBarTimeHists(self, hit):

        medians=self.muAna.GetMedianTime(hit, mode='globalcorrection') # Add mode option to functions for global correction
        if not medians: 
            print(f'Medians could not be computed for event: {self.M.EventNumber}, detID: {hit.GetDetectorID()}')
            return 
        averagetime=sum(medians.values())/len(medians)
        deltatime=medians['left'] - medians['right']
        
        detID=hit.GetDetectorID()
        subsystem, plane, bar= self.muAna.parseDetID(detID)
        averagebartimehistname=f'{10*subsystem+plane}_bar{bar}_averagetime'
        
        if not averagebartimehistname in self.hists:
            title=self.subsystemdict[subsystem]+' plane '+str(plane+1)+' bar '+str(bar+1)+' average of median '+self.state+' times from each side as a function of x_{predicted};x_{predicted} [cm];'+'#frac{1}{2}#times(t_{left}+t_{right}) [ns];Counts'
            self.hists[averagebartimehistname]=ROOT.TH2F(averagebartimehistname, title, 100, 0, 100, 200, -10, 10)
        self.hists[averagebartimehistname].Fill(self.xpred, averagetime)
        
        deltabartimehistname=f'{10*subsystem+plane}_bar{bar}_deltatime'
        if not deltabartimehistname in self.hists:
            title=self.subsystemdict[subsystem]+' plane '+str(plane+1)+' bar '+str(bar+1)+' difference of median '+self.state+' times from each side as a function of x_{predicted};x_{predicted} [cm];'+'t_{left} - t_{right} [ns];Counts'
            self.hists[deltabartimehistname]=ROOT.TH2F(deltabartimehistname, title, 100, 0, 100, 200, -10, 10)
        self.hists[deltabartimehistname].Fill(self.xpred, deltatime)

        for side in ('left', 'right'):
            averagebarsidetimehistname=f'{10*subsystem+plane}_bar{bar}_{side}averagetime'
            if not averagebarsidetimehistname in self.hists:
                title=self.subsystemdict[subsystem]+' plane '+str(plane+1)+' bar '+str(bar+1)+' average of median '+self.state+' times from each side as a function of x_{predicted};x_{predicted} [cm];'+'#frac{1}{2}#times(t_{left}+t_{right}) [ns];Counts'
                self.hists[averagebarsidetimehistname]=ROOT.TH2F(averagebarsidetimehistname, title, 100, 0, 100, 200, -10, 10)
            self.hists[averagebarsidetimehistname].Fill(self.xpred, medians[side])

    def WriteOutHistograms(self):
            
        testing_outfile=f'testing_{self.options.runNumber}_{self.options.nStart}_{self.options.nEvents}.root'
        testing_f=ROOT.TFile.Open(testing_outfile, 'recreate')
        
        additionalkeys=['averagetime', 'deltatime']
        for histname in self.M.h:
            if not any( [histname.find(additionalkey)>-1 for additionalkey in additionalkeys] ):
                
                hist=self.M.h[histname]
                testing_f.WriteObject(hist, hist.GetName(), 'kOverwrite')
            else:
                for additionalkey in additionalkeys:
                    if histname.find(additionalkey)>-1:
                        if not hasattr(testing_f, additionalkey): folder=testing_f.mkdir(additionalkey)
                        else: folder=testing_f.Get(additionalkey)
                        folder.cd()
                        hist=self.M.h[histname]
                        hist.Write(histname, 2) # The 2 means it will overwrite a hist of the same name
        testing_f.Close()
        print(f'{len(self.M.h)} histograms saved to {testing_outfile}')

    # def yresidual3(self, detID, pos, mom):
    def yresidual3(self, detID):
        self.MuFilter.GetPosition(detID,A,B)
        doca=self.Getyresidual(detID)
        s,p,b=self.muAna.parseDetID(detID)
        if s==3: return True
        key=10*s+p
        
        dy_min, dy_max = SystemAlignment.dycut(self.cutdists[f'dy_{key}'])
        if doca>dy_max or doca<dy_min: return False
        else: return True

    def Getyresidual(self, detID):
        s,p,b=self.muAna.parseDetID(detID)
        key=10*s+p
        if key>=30: return True
        # z=self.zPos['MuFilter'][key]

        pq = A-self.pos
        uCrossv= (B-A).Cross(self.mom)
        doca = pq.Dot(uCrossv)/uCrossv.Mag()    
        return doca

    def nSiPMscut(self, hit, nLeft, nRight):
        s,p,b=self.muAna.parseDetID(hit.GetDetectorID())
        if s==1: 
            if nLeft<6 or nRight<6: return False
        elif s==2:
            if nLeft<4 or nRight<4: return False
        elif s==3: 
            pass 
        return True
    
    def GetTimingDiscriminant(self):
        
        hits=self.M.eventTree.Digi_MuFilterHits
        US1hits=[h.GetDetectorID() for h in hits if all([h.GetDetectorID()//10000==2, self.muAna.parseDetID(h.GetDetectorID())[1]==0])]
        
        if len(US1hits)==0: return
        elif len(US1hits)==1:
             
            x=US1hits[0]
            tmp={h.GetDetectorID():h for h in hits}
            us1hit=tmp[x]
            
        elif len(US1hits)>1:
            docas={}
            for US1detID in US1hits:
                self.MuFilter.GetPosition(US1detID, A, B)
                docas[self.Getyresidual(US1detID)]=US1detID
            x=docas.pop(min(docas))
            tmp={h.GetDetectorID():h for h in hits}
            us1hit=tmp[x]  
        else: 
            print('shite')
            return

        averageUS1time=self.muAna.GetAverageTime(us1hit)
        if not averageUS1time: return -998

        DS3Haverage=self.muAna.GetDSHaverage(self.MuFilter, hits, mode='timingdiscriminant')

        return DS3Haverage-averageUS1time    

    def TimingDiscriminantCut(self, td):
        
        if self.timealignment=='old':
            if td<0: return False
            else: return True
        else: 
            if not 'timingdiscriminant' in self.cutdists:
                print(f'No timing discriminant histogram')
                return 
            hist=self.cutdists['timingdiscriminant']
            mean, stddev=hist.GetMean(), hist.GetStdDev()
            # print(f'td: {td}, mean +/- 2*std. dev.: {mean-2*stddev}, {mean+2*stddev}')
            if td < mean-2*stddev or td > mean+2*stddev: return False
            else: return True

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

    def TDS0cut(self,TDS0):
        if TDS0==-6237.5: return 0
        TDS0ns=TDS0*self.TDC2ns
        if TDS0ns<13.75 or TDS0ns>15.45: return 0
        return TDS0ns 
    
    def GetDistanceToSiPM(self,A,Ex):
        self.xpred=ROOT.TMath.Sqrt((A.x()-Ex.x())**2+(A.y()-Ex.y())**2+(A.z()-Ex.z())**2)
        # return ROOT.TMath.Sqrt((A.x()-Ex.x())**2+(A.y()-Ex.y())**2+(A.z()-Ex.z())**2)

    def muonvelocity(self, hits):
        
        # title='Average time recorded for single muon tracks in US system;z-position [cm];#bar{t_{bar}} [ns]' 
        finalhist=ROOT.TH1F('MuonVelocity',title,130,350,480)

        for idx,hit in enumerate(hits):
            SiPMtimes=[]
            detID=hit.GetDetectorID()
            s,p,b=self.muAna.parseDetID(detID)
            if s!=2:continue
            for ch in hit.GetAllTimes():
                SiPM,clock=ch
                SiPMtimes.append(t)
                      