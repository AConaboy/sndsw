import ROOT, os, csv, subprocess, cProfile, io, pstats
from argparse import ArgumentParser

from AnalysisFunctions import Analysis
from pathlib import Path
from array import array

deciles=[i/10 for i in range(11)]
A, B, locA, locB = ROOT.TVector3(), ROOT.TVector3(), ROOT.TVector3(), ROOT.TVector3()

ROOT.gStyle.SetTitleXOffset(0.7)
ROOT.gStyle.SetTitleYOffset(0.6)
for i in ('x', 'y', 'z','t'): ROOT.gStyle.SetTitleSize(0.04, i)

class SystemAlignment(ROOT.FairTask):

    def Init(self, options, monitor):
       
        self.state='corrected'
        options.state='corrected'
        self.runNr = str(options.runNumber).zfill(6)
        self.muAna=Analysis(options)
        self.timealignment=self.muAna.GetTimeAlignmentType(runNr=self.runNr)

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

        self.muAna.MakeTWCorrectionDict(self.TWCorrectionRun)
        self.muAna.MakeAlignmentParameterDict(self.TWCorrectionRun)
        self.muAna.Makecscintdict(self.TWCorrectionRun, 'corrected')
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

        self.freq=160.316E6
        self.TDC2ns=1E9/self.freq

        self.largeSiPMmap={0:0 ,1:1 ,3:2 ,4:3 ,6:4 ,7:5}
        self.verticalBarDict={0:1, 1:3, 2:5, 3:6}
        
        self.USxref=self.MuFilter.GetConfParF('MuFilter/UpstreamBarX')/2 # To be set as USbarlength/2.
        self.muAna.USxref=self.USxref
        self.Vetoxref=self.MuFilter.GetConfParF('MuFilter/VetoBarX')/2
        self.muAna.Vetoxref=self.Vetoxref
        self.muAna.xrefs={1:self.Vetoxref, 2:self.USxref}
        self.sides=('left', 'right')

        # self.correctionfunction=lambda ps, qdc: 1/sum( [ ps[i]*qdc**i for i in range(len(ps)) ] ) 
        # self.getaverage=lambda d, key, i:sum( [(d[key][i]) for k in range(2) ])/len(d)        

        self.hists=self.M.h
        # self.ChannelsWithNoParams=[]
        
        # self.CorrectionType=5
        # self.correctionparams = lambda ps : [y for x,y in enumerate(ps) if x%2==0]
        # self.correctionfunction = lambda ps, qdc : ps[3]*(qdc-ps[0])/( ps[1] + ps[2]*(qdc-ps[0])*(qdc-ps[0]) ) + ps[4]*(qdc-ps[0])

        self.gauss=ROOT.TF1('mygauss', 'abs([0])*abs([4])/(abs([2])*sqrt(2*pi))*exp(-0.5*((x-[1])/[2])**2)+abs([3])', 4)
        for idx,param in enumerate(('Constant', 'Mean', 'Sigma', 'y-offset', 'bin-width')): self.gauss.SetParName(idx, param)

        self.sigmatds0=0.263 # ns 

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
        self.TDS0, firedDSHbars=self.muAna.GetDSHaverage(self.MuFilter, event.Digi_MuFilterHits) # Now returns in ns
        if self.TDS0==-999.: print(f'Event {self.M.EventNumber} has no DS horizontal hits with 2 fired SiPMs')
        if not 'TDS0' in self.hists:
           self.hists['TDS0']=ROOT.TH1F('TDS0','Average time of DS horizontal bars;DS horizontal average time [ns];Counts', 200, 0, 50)
        self.hists['TDS0'].Fill(self.TDS0)            

        # if self.options.timingdiscriminantcut:
        td = self.GetTimingDiscriminant() # Require that US1 TDC average is less than the DSH TDC average to ensure forward travelling track
        if not td: return
        if not self.TimingDiscriminantCut(td): return
        
        ### Slope cut
        if not SystemAlignment.slopecut(self.mom): return

        for hit in event.Digi_MuFilterHits:

            nLeft, nRight=self.muAna.GetnFiredSiPMs(hit)
            
            detID=hit.GetDetectorID()
            s,p,b=self.muAna.parseDetID(detID)
            
            if s==3: continue

            # Only investigate veto hits if there is a Scifi track!
            if not inVeto and s==1: continue

            # This function also calls MuFilter.GetPosition with the detID so it's important even for the DS, where the 
            # y-residual cut is not applied.
            self.MuFilter.GetPosition(detID,A,B)
            self.MuFilter.GetLocalPosition(detID, locA, locB)
            if not self.yresidual3(detID): continue

            zEx=self.zPos['MuFilter'][s*10+p]
            lam=(zEx-self.pos.z())/self.mom.z()
            self.Ex=ROOT.TVector3(self.pos.x()+lam*self.mom.x(), self.pos.y()+lam*self.mom.y(), self.pos.z()+lam*self.mom.z())

            ### Particle ToF for old time alignment data
            if self.timealignment=='old': self.particleToFcorrection=self.ParticleToFcorrection()

            self.GetDistanceToSiPM(A)

            channels_t=hit.GetAllTimes()
            channels_qdc=hit.GetAllSignals()

            if s!=3:
                self.FillSiPMHists(hit)
                # self.FillBarTimeHists(hit)

                if self.options.CrossTalk:
                    self.XTHists(hit)

    def FillSiPMHists(self, hit):

        detID=hit.GetDetectorID()

        # qdcs, clocks = hit.GetAllSignals(), hit.GetAllTimes()

        alignedtimes=self.muAna.GetAlignedTimes(hit, self.pred, mode='unaligned')

        for ch in alignedtimes:
            SiPM, time=ch
            fixed_ch=f'{detID}_{SiPM}'

            ReadableDetID=self.muAna.MakeHumanReadableFixedCh(fixed_ch)
            SiPMtime=f'AlignedSiPMtime_{fixed_ch}'
            if not SiPMtime in self.hists:
                title=f'DS horizontal average relative, TW + ToF corrected + DS aligned SiPM time'
                splittitle='#splitline{'+ReadableDetID+'}{'+title+'}'
                axestitles='t^{0}_{DS} - t_{SiPM '+str(SiPM)+'}^{tw corr + aligned} [ns];Counts'
                fulltitle=splittitle+';'+axestitles
                if self.timealignment=='old': self.hists[SiPMtime]=ROOT.TH1F(SiPMtime,fulltitle, 400, 0, 20)
                else: self.hists[SiPMtime]=ROOT.TH1F(SiPMtime,fulltitle, 400, -10, 10)

            alignmentparameter=self.muAna.alignmentparameters[fixed_ch]
            print(f'{fixed_ch}, {self.TDS0 - time - alignmentparameter[0] }')
            # alignedtime=time - alignmentparameter[0]
            self.hists[SiPMtime].Fill(self.TDS0 - time - alignmentparameter[0])

    def FillBarTimeHists(self, hit):

        """
        I will keep the Analysis class methods to return the 
        SiPM times relative to the event-t0. In the methods in this SystemAlignment class
        I will shift them to be relative to tds0.
        """
        
        ### medians returns a dictionary: {'left': float, 'right':float}
        # medians=self.muAna.GetMedianTime(hit, mode='alignment')
        alignedtimes=self.muAna.GetAlignedTimes(hit, self.pred)
        
        if medians==-999: return
        if not medians: 
            # print(f'Medians could not be computed for event: {self.M.EventNumber}, detID: {hit.GetDetectorID()}')
            return 

        # averagetime=sum(medians.values())/len(medians)
        averagetime=1/2 * sum([self.TDS0-i for i in medians.values()])

        deltatime=self.TDS0-medians['left'] - self.TDS0-medians['right']
        
        detID=hit.GetDetectorID()
        subsystem, plane, bar= self.muAna.parseDetID(detID)
        
        averagebartimehistname=f'{10*subsystem+plane}_bar{bar}_averagetime'
        if not averagebartimehistname in self.hists:
            title=self.subsystemdict[subsystem]+' plane '+str(plane+1)+' bar '+str(bar+1)+' average of median '+self.state+' times from each side as a function of x_{predicted};x_{predicted} [cm];'+'#frac{1}{2}#times(t^{tw corr}_{left}+t^{tw corr}_{right}) [ns];Counts'
            self.hists[averagebartimehistname]=ROOT.TH2F(averagebartimehistname, title, 100, 0, 100, 200, -10, 10)
        self.hists[averagebartimehistname].Fill(self.pred, self.TDS0-averagetime)
        
        deltabartimehistname=f'{10*subsystem+plane}_bar{bar}_deltatime'
        if not deltabartimehistname in self.hists:
            title=self.subsystemdict[subsystem]+' plane '+str(plane+1)+' bar '+str(bar+1)+' difference of median '+self.state+' times from each side as a function of x_{predicted};x_{predicted} [cm];'+'t^{tw corr}_{left} - t^{tw corr}_{right} [ns];Counts'
            self.hists[deltabartimehistname]=ROOT.TH2F(deltabartimehistname, title, 100, 0, 100, 200, -10, 10)
        self.hists[deltabartimehistname].Fill(self.pred, deltatime)

        for side in ('left', 'right'):
            averagebarsidetimehistname=f'{10*subsystem+plane}_bar{bar}_{side}sidetime'
            if not averagebarsidetimehistname in self.hists:
                title=self.subsystemdict[subsystem]+' plane '+str(plane+1)+' bar '+str(bar+1)+' median '+self.state+' times from '+side+' side as a function of x_{predicted};x_{predicted} [cm];'+'Median of t^{tw corr}_{'+side+'} [ns];Counts'
                self.hists[averagebarsidetimehistname]=ROOT.TH2F(averagebarsidetimehistname, title, 100, 0, 100, 200, -10, 10)
            self.hists[averagebarsidetimehistname].Fill(self.pred, self.TDS0- medians[side])

    def XTHists(self, hit):

        clocks, qdcs =hit.GetAllTimes(), hit.GetAllSignals()
        times={'left':{}, 'right':{}}
        detID=hit.GetDetectorID()
        s,p,b = self.muAna.parseDetID(detID)

        for i in clocks: 
            SiPM, clock=i
            fixed_ch=self.muAna.MakeFixedCh((s,p,b,SiPM))
            qdc=self.muAna.GetChannelVal(SiPM, qdcs)
            correctedtime=self.muAna.MuFilterCorrectedTime(fixed_ch, qdc, clock, self.pred)
            if not correctedtime: continue
            side=self.muAna.GetSide(f'{detID}_{SiPM}')
            times[side][SiPM]=correctedtime

        for side in times:
            for i in times[side]:
                for j in times[side]:
                    if i==j: continue
                    permutation_id = f'{i}{j}'
                    xthistname=f'timingxt_{detID}_SiPMs{i}-{j}'
                    if not xthistname in self.hists:
                        title=f'Timing correlation, {self.subsystemdict[s]} plane {p+1} bar {b+1}'
                        subtitle=f'{side} SiPMs {i+1} and {j+1}'
                        splittitle='#splitline{'+title+'}{'+subtitle+'}'
                        axestitles='t_{0}^{DS} - t_{SiPM '+str(i)+'}^{tw corr + aligned} [ns];t_{0}^{DS} - t_{SiPM '+str(j)+'}^{tw corr + aligned} [ns]'
                        fulltitle=splittitle+';'+axestitles
                        if self.timealignment=='old': self.hists[xthistname]=ROOT.TH1F(xthistname,fulltitle, 100, 0, 20, 100, 0, 20)
                        else: self.hists[xthistname]=ROOT.TH2F(xthistname,fulltitle, 100, -5, 15, 100, -5, 15)
                    ti, tj= self.TDS0-times[side][i],self.TDS0-times[side][j]
                    self.hists[xthistname].Fill(ti, tj)

    def WriteOutHistograms(self):

        outfilename=f'{self.outpath}splitfiles/run{self.runNr}/SystemAlignment/SystemAlignment_{self.options.nStart}.root'
        if os.path.exists(outfilename): outfile=ROOT.TFile.Open(outfilename, 'recreate')
        else: outfile=ROOT.TFile.Open(outfilename, 'create')

        additionalkeys=['averagetime', 'deltatime', 'sidetime', 'AlignedSiPMtime', 'timingxt']
        for h in self.M.h:
            if h=='TDS0':
                outfile.WriteObject(self.M.h[h], self.M.h[h].GetName(),'kOverwrite')

            if len(h.split('_'))==3:
                for additionalkey in additionalkeys:
                    if h.find(additionalkey)==-1: continue
                    planekey, bar, key = h.split('_')
                    
                    hist=self.M.h[h]
                    if not hasattr(outfile, additionalkey): folder=outfile.mkdir(additionalkey)
                    else: folder=outfile.Get(additionalkey)
                    folder.cd()
                    hist.Write(h, 2) # The 2 means it will overwrite a hist of the same name                

        outfile.Close()
        print(f'{len(self.M.h)} histograms saved to {outfilename}')        

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
    
    def GetDistanceToSiPM(self,A):
        self.pred=ROOT.TMath.Sqrt((A.x()-self.Ex.x())**2+(A.y()-self.Ex.y())**2+(A.z()-self.Ex.z())**2)

    def ParticleToFcorrection(self, start=('Scifi','10')):
        zPos_start=self.zPos[start[0]][int(start[1])]
        lam=(zPos_start-self.pos.z())/self.mom.z()
        Ex_scifi=ROOT.TVector3(self.pos.x()+lam*self.mom.x(), self.pos.y()+lam*self.mom.y(), self.pos.z()+lam*self.mom.z())
        R_sq = (self.Ex.x()-Ex_scifi.x())**2 + (self.Ex.y()-Ex_scifi.y())**2 + (self.Ex.z()-Ex_scifi.z())**2
        R=ROOT.TMath.Sqrt(R_sq)

        tofcorrection=R/30
        return tofcorrection         

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
                      