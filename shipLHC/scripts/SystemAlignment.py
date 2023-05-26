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
        statedict={'zeroth':'uncorrected', 'ToF':'uncorrected', 'TW':'corrected'}
        self.state=statedict[options.mode]

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
        self.correctionparams=lambda ps : [y for x,y in enumerate(ps) if x%2==0]
        self.correctionfunction = lambda ps, qdc : ps[3]*(qdc-ps[0])/( ps[1] + ps[2]*(qdc-ps[0])*(qdc-ps[0]) ) 

        if options.TWCorrectionRun!=None: self.TWCorrectionRun=str(options.TWCorrectionRun).zfill(6)
        else: self.TWCorrectionRun=self.runNr

    def GetEntries(self):
        return self.eventTree.GetEntries()

    def ExecuteEvent(self, event):

        tracks={}
        Reco_MuonTracks=self.M.Reco_MuonTracks
        inVeto, inDS=False, False
        for i,track in enumerate(Reco_MuonTracks):
            if track.GetUniqueID()==1 and track.getFitStatus().isFitConverged(): 
                inVeto=True
                tracks[1]=Reco_MuonTracks[i]
            if track.GetUniqueID()==3: 
                # if not track.getFitStatus().isFitConverged(): return
                inDS=True
                tracks[2]=Reco_MuonTracks[i]
                DStrack=tracks[2]
        if not inDS: return

        fstate=DStrack.getFittedState()
        fitStatus=DStrack.getFitStatus()
        if not fitStatus.isFitConverged(): return
        # fstates=[i:tracks[i].getFittedState() for i in tracks]
        self.pos=fstate.getPos()
        # posvectors=[i:fstates[i].getPos() for i in fstates]
        self.mom=fstate.getMom()
        # momvectors=[i:fstates[i].getMom() for i in fstates]
        self.trackchi2NDF=fitStatus.getChi2()/fitStatus.getNdf()+1E-10

        # Get DS event t0
        TDS0=self.muAna.GetDSHaverage(self.MuFilter, event.Digi_MuFilterHits) # Now returns in ns
        if TDS0==-999.: print(f'Event {self.M.EventNumber} has no DS horizontal hits with 2 fired SiPMs')
        if not 'TDS0' in self.hists:
           self.hists['TDS0']=ROOT.TH1F('TDS0','Average time of DS horizontal bars;DS horizontal average time [ns];Counts', 200, 0, 50)
        self.hists['TDS0'].Fill(TDS0)

        # if self.options.timingdiscriminantcut:
        td = self.GetTimingDiscriminant() # Require that US1 TDC average is less than the DSH TDC average to ensure forward travelling track
        if not self.TimingDiscriminantCut(td): return

        ### Slope cut
        if not TimeWalk.slopecut(self.mom): return

        for hit in event.Digi_MuFilterHits:
            nLeft, nRight=self.muAna.GetnFiredSiPMs(hit)
            
            # if not self.nSiPMscut(hit, nLeft, nRight): continue

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

            self.FillHists(self, hit)

            # for channel in channels_t:
            #     SiPM,clock=channel
            #     qdc=self.muAna.GetChannelVal(SiPM, channels_qdc)
            #     if qdc==-999.: continue
            #     fixed_ch=f'{detID}_{SiPM}'

            #     self.polyparams=self.muAna.GetPolyParams(self.TWCorrectionRun, fixed_ch)
            #     if self.polyparams==-999.:
            #         if fixed_ch not in self.ChannelsWithNoParams: self.ChannelsWithNoParams.append(fixed_ch)


    def FillHists(self, hit):

        medians=self.muAna.GetMedianTime(hit, state='globalcorrection') # Add mode option to functions for global correction
        if medians==-999: 
            print(f'Medians could not be computed for {hit.GetDetectorID()}')
            return 
        averagetime=sum(medians)/len(medians)

        averagebartimehistname=f'{10*subsystem+plane}_bar{bar}_averagetime'
        if not averagebartimehistname in self.hists:
            title=self.subsystemdict[subsystem]+' plane '+str(plane+1)+' bar '+str(bar+1)+' average of median '+self.state+' time from each side as a function of x_{predicted};x_{predicted} [cm];'+'#frac{1}{2}#times(t_{left}+t_{right}) [ns];Counts'
            self.hists[averagebartimehistname]=ROOT.TH2F(averagebartimehistname, title, 100,0,100,250,0,25)
        self.hists[averagebartimehistname].Fill(self.xpred, )

        


    def WriteOutHistograms(self):
        testing_outfile=f'testing_{self.options.runNumber}_{self.options.nStart}_{self.options.nEvents}.root'
        testing_f=ROOT.TFile.Open(testing_outfile, 'update')                    

        for h in self.M.h:
            if self.options.debug:
                if h.find('channelhitrate')!=-1:
                    hist=self.M.h[h]
                    testing_f.WriteObject(hist, f'TimeWalk-{h}')
            if len(h.split('_'))==4:
                histkey,detID,SiPM,state=h.split('_')
                fixed_ch='_'.join((detID,SiPM))
                hist=self.M.h[h]
                outpath=f'{self.outpath}/splitfiles/run{self.runNr}/{fixed_ch}/'
                path_obj=Path(outpath)
                path_obj.mkdir(parents=True, exist_ok=True)
                
                outfile=f'timewalk_{fixed_ch}_{self.options.nStart}.root'
                f=ROOT.TFile.Open(outpath+outfile, 'update')
                f.WriteObject(hist, hist.GetName(), 'kOverwrite')
                if self.options.mode=='zeroth':
                    f.WriteObject(self.hists['TDS0'], 'TDS0', 'kOverwrite')
                f.Close()
        print(f'{len(self.M.h)} histograms saved to {self.outpath}splitfiles/run{self.runNr}/fixed_ch')
        
        if self.options.debug:
            print(f'Channel hit rate hists written to {testing_outfile}')
            testing_f.Close()

    # def yresidual3(self, detID, pos, mom):
    def yresidual3(self, detID):
        self.MuFilter.GetPosition(detID,A,B)
        doca=self.Getyresidual(detID)
        s,p,b=self.muAna.parseDetID(detID)
        if s==3: return True
        key=10*s+p
        
        dy_min, dy_max = TimeWalk.dycut(self.cutdists[f'dy_{key}'])
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
        
        if len(US1hits)==0: return -999
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
            return -999

        averageUS1time=self.muAna.GetAverageTime(us1hit)
        if not averageUS1time: return -998

        DS3Haverage=self.muAna.GetDSHaverage(self.MuFilter, hits, mode='timingdiscriminant')

        return DS3Haverage-averageUS1time    

    def TimingDiscriminantCut(self, td):
        timingalignment=self.muAna.GetTimeAlignmentType(runNr=self.runNr)
        if timingalignment=='old':
            if td<0: return False 
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
        if TDS0==-999. or TDS0==-998. or TDS0==-6237.5: return 0
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
                      