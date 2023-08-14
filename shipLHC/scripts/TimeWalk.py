import ROOT, os, csv, subprocess, cProfile, io, pstats, json
from argparse import ArgumentParser
from AnalysisFunctions import Analysis
from pathlib import Path
from array import array

deciles=[i/10 for i in range(11)]

ROOT.gStyle.SetTitleXOffset(0.7)
ROOT.gStyle.SetTitleYOffset(0.6)
for i in ('x', 'y', 'z','t'): ROOT.gStyle.SetTitleSize(0.04, i)

class TimeWalk(ROOT.FairTask):

    def Init(self, options, monitor):
        
        statedict={'zeroth':'uncorrected', 'ToF':'uncorrected',
                    'TW':'corrected', 'res':'corrected', 
                    'systemalignment':'corrected'}
        self.state=statedict[options.mode]
        
        self.A, self.B = ROOT.TVector3(), ROOT.TVector3()
        self.locA, self.locB =ROOT.TVector3(), ROOT.TVector3()
        
        self.muAna=Analysis(options)
        self.M=monitor
        self.options=options
        if self.options.path.find('commissioning/TI18')>0:
            self.outpath=options.afswork+'-commissioning/'
            self.path='TI18'
        elif self.options.path.find('physics/202')>0: 
            self.outpath=options.afswork+'-physics2022/'
            self.path='TI18'
        elif self.options.path.find('H8')>0: 
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
        self.muAna.zPos=self.zPos

        self.timealignment=self.muAna.GetTimeAlignmentType(runNr=self.runNr)

        self.freq=160.316E6
        self.TDC2ns=1E9/self.freq

        self.largeSiPMmap={0:0 ,1:1 ,3:2 ,4:3 ,6:4 ,7:5}
        self.verticalBarDict={0:1, 1:3, 2:5, 3:6}
        
        self.USxref=self.MuFilter.GetConfParF('MuFilter/UpstreamBarX')/2 # To be set as USbarlength/2.
        self.Vetoxref=self.MuFilter.GetConfParF('MuFilter/VetoBarX')/2
        self.xrefs={1:self.Vetoxref, 2:self.USxref}
        self.muAna.xrefs=self.xrefs

        if options.debug: self.trackevents=[]
        
        # self.correctionfunction=lambda ps, qdc: 1/sum( [ ps[i]*qdc**i for i in range(len(ps)) ] ) 
        self.getaverage=lambda d, key, i:sum( [(d[key][i]) for k in range(2) ])/len(d)        

        self.hists=self.M.h

        if self.options.debug:
            systemhistcriteria={None:None, '1DSb/p':'1 DS bar per plane', '1USb/p':'1 US bar per plane'}
            for s in (1, 2):
                for criterion in systemhistcriteria.keys():
                    if not criterion:
                        title=self.subsystemdict[s]+' channel hit rate;Channel;Counts'
                        self.hists[f'channelhitrate-{self.subsystemdict[s]}']=ROOT.TH1F(f'channelhitrate-{self.subsystemdict[s]}', title, self.nchs[s]+1, 0, self.nchs[s])
                        
                        title='#splitline{'+self.subsystemdict[s]+' channel hit rate}{with timing discriminant cut};Channel;Counts'
                        self.hists[f'channelhitrate-{self.subsystemdict[s]}-tdcut']=ROOT.TH1F(f'channelhitrate-{self.subsystemdict[s]}-tdcut', title, self.nchs[s]+1, 0, self.nchs[s])

                        title='#splitline{'+self.subsystemdict[s]+' channel hit rate}{with dy cut};Channel;Counts'
                        self.hists[f'channelhitrate-{self.subsystemdict[s]}-dycut']=ROOT.TH1F(f'channelhitrate-{self.subsystemdict[s]}-dycut', title, self.nchs[s]+1, 0, self.nchs[s])
                    else: 
                        
                        title='#splitline{'+self.subsystemdict[s]+' channel hit rate}{with '+systemhistcriteria[criterion]+'};Channel;Counts'
                        self.hists[f'channelhitrate-{self.subsystemdict[s]}-{criterion}']=ROOT.TH1F(f'channelhitrate-{self.subsystemdict[s]}-{criterion}', title, self.nchs[s]+1, 0, self.nchs[s])
        
        # if options.mode in ('TW', 'res'):
        if not options.CorrectionType: self.CorrectionType=4
        else: self.CorrectionType=options.CorrectionType

        self.correctionparams=lambda ps : [y for x,y in enumerate(ps) if x%2==0]
        
        if options.CorrectionType==1: self.correctionfunction = lambda ps, qdc : 1/sum( [ ps[i]*qdc**i for i in range(len(ps)) ] ) 
        elif options.CorrectionType==4: self.correctionfunction = lambda ps, qdc : ps[3]*(qdc-ps[0])/( ps[1] + ps[2]*(qdc-ps[0])*(qdc-ps[0]) )
        elif options.CorrectionType==5: self.correctionfunction = lambda ps, qdc : ps[3]*(qdc-ps[0])/( ps[1] + ps[2]*(qdc-ps[0])*(qdc-ps[0]) ) + ps[4]*(qdc-ps[0])

        ### If no time-walk correction run is provided. Set the default correction run depending on time alignment of the data set
        if options.TWCorrectionRun==None:
            if self.timealignment=='old': self.TWCorrectionRun=str(5097).zfill(6)
            elif self.timealignment=='new': self.TWCorrectionRun=str(5408).zfill(6)
            elif self.timealignment=='new+LHCsynch': self.TWCorrectionRun=str(5999).zfill(6)
        ### Check that the time-walk correction run selected has the same time alignment has the data set
        else:
            if self.muAna.GetTimeAlignmentType(runNr=options.TWCorrectionRun) != self.timealignment:
                if self.timealignment=='old': self.TWCorrectionRun=str(5097).zfill(6)
                elif self.timealignment=='new': self.TWCorrectionRun=str(5408).zfill(6)
                elif self.timealignment=='new+LHCsynch': self.TWCorrectionRun=str(5999).zfill(6)
            else: self.TWCorrectionRun=str(options.TWCorrectionRun).zfill(6)
        self.muAna.TWCorrectionRun=self.TWCorrectionRun

        self.cutdists=self.muAna.GetCutDistributions(self.TWCorrectionRun, ('dy', 'timingdiscriminant'))
        ### Incase selection criteria distributions not made for the TWCorrectionRun
        if not self.cutdists:
            self.TWCorrectionRun=str(5408).zfill(6)
            self.muAna.TWCorrectionRun=self.TWCorrectionRun
            self.cutdists=self.muAna.GetCutDistributions(self.TWCorrectionRun, ('dy', 'timingdiscriminant'))

        ### Dictionary of cscint values using time-walk correction run
        if options.mode != 'zeroth':
            self.muAna.Makecscintdict(self.TWCorrectionRun, state=self.state)

        if options.mode in ('TW', 'res', 'systemalignment'):
            self.muAna.MakeTWCorrectionDict(self.TWCorrectionRun)

        if options.mode == 'systemalignment':
            self.muAna.MakeAlignmentParameterDict(self.TWCorrectionRun)
            from systemalignment import SystemAlignment
            self.sa=SystemAlignment(options, self)
        
        self.scifi1z=self.zPos['Scifi'][10]

    def GetEntries(self):
        return self.eventTree.GetEntries()

    def ExecuteEvent(self, event):

        if not hasattr(self.muAna, 'task'): self.muAna.SetTask(self)

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

        ### If there are more than 1 DS track, take the track with the lowest chi2/Ndf
        if len(tracks[3])==1: dstrack= tracks[3][0]
        else: 
            tmp={i:i.getFitStatus().getChi2()/i.getFitStatus().getNdf() for i in tracks[3]}
            dstrack=tmp[ min(tmp) ]
        
        fitStatus=dstrack.getFitStatus()
        fstate=dstrack.getFittedState()

        self.pos=fstate.getPos()
        self.mom=fstate.getMom()
        self.trackchi2NDF=fitStatus.getChi2()/fitStatus.getNdf()+1E-10

        # if self.options.debug:
        #     self.FillChannelRateHists()

        # Get DS event t0
        x=self.muAna.GetDSHaverage(event.Digi_MuFilterHits) # Now returns in ns
        if x==-999.: 
            print(f'Event {self.M.EventNumber} has no DS horizontal hits with 2 fired SiPMs')
            return
        self.TDS0, firedDSHbars=x
        if not 'TDS0' in self.hists:
           self.hists['TDS0']=ROOT.TH1F('TDS0','Average time of DS horizontal bars;DS horizontal average time [ns];Counts', 200, 0, 50)
        self.hists['TDS0'].Fill(self.TDS0)

        ### Timing discriminant cut
        td = self.muAna.GetTimingDiscriminant() # Require that US1 TDC average is less than the DSH TDC average to ensure forward travelling track
        if not self.TimingDiscriminantCut(td): return
        
        ### Slope cut
        if not self.slopecut(): return
    
        for hit in event.Digi_MuFilterHits:
            nLeft, nRight=self.muAna.GetnFiredSiPMs(hit)
            
            # if not self.nSiPMscut(hit, nLeft, nRight): continue

            detID=hit.GetDetectorID()
            s,p,b=self.muAna.parseDetID(detID)

            # Only investigate veto hits if there is a Scifi track!
            if not inVeto and s==1: continue

            # This function also calls MuFilter.GetPosition with the detID so it's important even for the DS, where the 
            # y-residual cut is not applied.
            self.MuFilter.GetPosition(detID,self.A,self.B)
            self.MuFilter.GetLocalPosition(detID, self.locA, self.locB)
            if not self.yresidual3(detID): continue

            zEx=self.zPos['MuFilter'][s*10+p]
            lam=(zEx-self.pos.z())/self.mom.z()
            Ex=ROOT.TVector3(self.pos.x()+lam*self.mom.x(), self.pos.y()+lam*self.mom.y(), self.pos.z()+lam*self.mom.z())

            self.pred=self.GetDistanceToSiPM(Ex)

            if self.options.debug: self.trackevents.append(self.M.EventNumber)

            channels_t=hit.GetAllTimes()
            channels_qdc=hit.GetAllSignals()

            if self.options.mode=='systemalignment':
                self.sa.FillSiPMHists(hit)
                continue

            for channel in channels_t:
                SiPM,clock=channel
                qdc=self.muAna.GetChannelVal(SiPM, channels_qdc)
                if qdc==-999.: continue
                fixed_ch=f'{detID}_{SiPM}'
                if self.options.mode=='zeroth': self.zeroth(fixed_ch, clock, qdc)
                elif self.options.mode=='ToF': self.ToF(fixed_ch, clock, qdc)
                elif self.options.mode=='TW': self.TW(fixed_ch, clock, qdc, meantimecorrection=False)
                elif self.options.mode=='res': self.res(fixed_ch, clock, qdc)

    def FillChannelRateHists(self):
        hits=self.M.eventTree.Digi_MuFilterHits
        OneHitPerUS=self.muAna.OneHitPerSystem(hits, (2,))
        OneHitPerDS=self.muAna.OneHitPerSystem(hits, (3,))
        for hit in hits:
            detID=hit.GetDetectorID()
            s,p,b=self.muAna.parseDetID(detID)

            if s==3: continue
            for x in hit.GetAllTimes():
                SiPM, clock = x
                channel=self.muAna.GetSiPMNumberInSystem_PCBbyPCB(detID, SiPM)
                if OneHitPerUS: self.hists[f'channelhitrate-{self.subsystemdict[s]}-1USb/p'].Fill(channel)
                if OneHitPerDS: self.hists[f'channelhitrate-{self.subsystemdict[s]}-1DSb/p'].Fill(channel)
                self.hists[f'channelhitrate-{self.subsystemdict[s]}'].Fill(channel)
                td=self.GetTimingDiscriminant()
                if self.TimingDiscriminantCut(td): self.hists[f'channelhitrate-{self.subsystemdict[s]}-tdcut'].Fill(channel)
                if self.yresidual3(detID): self.hists[f'channelhitrate-{self.subsystemdict[s]}-dycut'].Fill(channel)

    def zeroth(self,fixed_ch,clock,qdc):
        hists=self.hists
        
        detID, SiPM = ( int(fixed_ch.split('_')[i]) for i in range(2) )
        s,p,b=self.muAna.parseDetID(detID)

        correctedtime=clock*self.TDC2ns
        t_rel=self.TDS0-correctedtime
        if not self.muAna.DSVcheck(detID):  dtvpred=f'dtvxpred_{fixed_ch}_{self.state}'
        else: dtvpred=f'dtvypred_{fixed_ch}_{self.state}'
        
        ReadableFixedCh=self.muAna.MakeHumanReadableFixedCh(fixed_ch)
    
        if not dtvpred in hists:
            coord='x' if not self.muAna.DSVcheck(detID) else 'y'
            title='{No time-walk correction t_{0}^{DS}-t^{uncorr}_{SiPM} v '+coord+'-position}'
            splittitle='#splitline{'+ReadableFixedCh+'}'+title
            axestitles=coord+'_{predicted} [cm];t_{0}^{DS}-t^{uncorr}_{SiPM} [ns]'
            fulltitle=splittitle+';'+axestitles
            hists[dtvpred]=ROOT.TH2F(dtvpred,fulltitle,110,-10,100, 800, -20, 20.)

        SiPMtime=f'tSiPM_{fixed_ch}_{self.state}'
        if not SiPMtime in hists:
            title='{No time-walk correction t_{SiPM}}'
            splittitle='#splitline{'+ReadableFixedCh+'}'+title
            axestitles='t_{SiPM} [ns];Counts'
            fulltitle=splittitle+';'+axestitles
            if self.timealignment=='old': hists[SiPMtime]=ROOT.TH1F(SiPMtime,fulltitle, 400, 0, 20)
            else: hists[SiPMtime]=ROOT.TH1F(SiPMtime,fulltitle, 400, -10, 10)
        
        attlen=f'attlen_{fixed_ch}_{self.state}'
        if not attlen in hists:
            coord='x' if not self.muAna.DSVcheck(detID) else 'y'
            title='{Predicted position against QDC_{SiPM}}'
            splittitle='#splitline{'+ReadableFixedCh+'}'+title
            axestitles=coord+'_{predicted} [cm];QDC_{SiPM} [a.u]'
            fulltitle=splittitle+';'+axestitles 
            hists[attlen]=ROOT.TH2F(attlen,fulltitle, 110, -10, 110, 200, 0., 200)

        if self.options.debug and s==3 and self.muAna.DSVcheck(detID):
            verticalplanehitrate=f'verticalplanehitrate_plane{p}'
            if not verticalplanehitrate in hists:
                title=f'Vertical plane hit rate for plane {p+1};Bar number;Counts'
                hists[verticalplanehitrate]=ROOT.TH1F(verticalplanehitrate,title, 120, 0, 119)
            hists[verticalplanehitrate].Fill(b)
        
        self.hists[dtvpred].Fill(self.pred,t_rel)
        self.hists[SiPMtime].Fill(correctedtime)
        self.hists[attlen].Fill(self.pred, qdc)

    def ToF(self, fixed_ch, clock, qdc):        
        hists=self.hists
        ReadableFixedCh=self.muAna.MakeHumanReadableFixedCh(fixed_ch)

        # cdata=self.muAna.Getcscint(self.runNr, fixed_ch, self.state)
        if not fixed_ch in self.muAna.cscintvalues:return
        # cdata=self.muAna.cscintvalues[fixed_ch]
        cdata=self.systemobservables[fixed_ch]['cscint']['uncorrected']
    
        s, SiPM=int(fixed_ch[0]), int(fixed_ch.split('_')[-1])
        correctedtime=self.muAna.correct_ToF(SiPM, clock, self.pred, cdata, self.xrefs[s])[1]
        dtvqdc=f'dtvqdc_{fixed_ch}_{self.state}'
        if not dtvqdc in hists:
            subtitle='{No time-walk correction t_{0}^{DS}-t^{uncorr}_{SiPM} v QDC_{SiPM}};QDC_{SiPM} [a.u];t_{0}^{DS}-t^{uncorr}_{SiPM} [ns]'
            title='#splitline{'+ReadableFixedCh+'}'+subtitle
            hists[dtvqdc]=ROOT.TH2F(dtvqdc, title, 200, 0, 200, 800, -20, 20)
        t_rel=self.TDS0-correctedtime
        self.hists[dtvqdc].Fill(qdc,t_rel)

    def TW(self, fixed_ch, clock, qdc, meantimecorrection=False):
        hists=self.hists
        ReadableFixedCh=self.muAna.MakeHumanReadableFixedCh(fixed_ch)
        
        s=int(fixed_ch[0])
        detID, SiPM=int(fixed_ch.split('_')[0]), int(fixed_ch.split('_')[1])
        time=clock*self.TDC2ns
        
        ### Check if fixed_ch has tw corr params and a cscint value
        if any([fixed_ch not in self.muAna.twparameters, fixed_ch not in self.muAna.cscintvalues]): return 

        cdata=self.muAna.cscintvalues[fixed_ch]
        twparams=self.muAna.twparameters[fixed_ch]

        ToFcorrectedtime=self.muAna.correct_ToF(SiPM, clock, cdata)[1]
        twcorrection = self.correctionfunction(twparams, qdc)

        ### TW corrected time then ToF & TW corrected time
        TWcorrectedtime=time+twcorrection
        ToFTWcorrectedtime=ToFcorrectedtime+twcorrection

        ### Times wrt to DS horizontal average
        ToFTWt_rel = self.TDS0-ToFTWcorrectedtime
        TWt_rel = self.TDS0-TWcorrectedtime

        ### Make histograms
        if not self.muAna.DSVcheck(detID):  dtvxpred=f'dtvxpred_{fixed_ch}_{self.state}'
        else: dtvxpred=f'dtvypred_{fixed_ch}_{self.state}'
        if not dtvxpred in hists:
            coord='x' if not self.muAna.DSVcheck(detID) else 'y'
            subtitle='{Time-walk corrected t_{0}^{DS}-t^{tw corr}_{SiPM} v '+coord+'-position w/ run'+str(self.TWCorrectionRun)+'};'+coord+'_{predicted} [cm];t_{0}^{DS}-t^{tw corr}_{SiPM} [ns]'
            title='#splitline{'+MakeHumanReadableFixedCh+'}'+subtitle
            hists[dtvxpred]=ROOT.TH2F(dtvxpred,title,110,-10,100, 800, -20, 20.)

        SiPMtime=f'tSiPM_{fixed_ch}_corrected'
        if not SiPMtime in hists:
            subtitle='{Time-walk corrected SiPM time w/ run'+str(self.TWCorrectionRun)+'};t^{tw corr}_{SiPM '+str(SiPM)+'} [ns];Counts'
            title='#splitline{'+MakeHumanReadableFixedCh+'}'+subtitle
            # hists[SiPMtime]=ROOT.TH1F(SiPMtime,title, 400, -10, 10)
            if self.timealignment=='old': hists[SiPMtime]=ROOT.TH1F(SiPMtime,title, 400, 0, 20)
            else: hists[SiPMtime]=ROOT.TH1F(SiPMtime,title, 400, -10, 10) 

        ### Fill histograms 
        hists[dtvxpred].Fill(self.pred, TWt_rel)
        hists[SiPMtime].Fill(TWcorrectedtime)
        
    def res(self, fixed_ch, clock, qdc):
        
        hists=self.hists
        ReadableFixedCh=self.muAna.MakeHumanReadableFixedCh(fixed_ch)
        detID, SiPM = int(fixed_ch.split('_')[0]), int(fixed_ch.split('_')[1])
        ### Check if fixed_ch has tw corr params and a cscint value
        if any([fixed_ch not in self.muAna.twparameters, fixed_ch not in self.muAna.cscintvalues]): return 

        cdata=self.muAna.cscintvalues[fixed_ch]
        twparams=self.muAna.twparameters[fixed_ch]

        ToFtime=self.muAna.correct_ToF(fixed_ch, clock, self.pred)[1]
        twcorrection = self.correctionfunction(twparams, qdc)

        ### TW corrected time then ToF & TW corrected time
        ToFTWtime=ToFtime+twcorrection

        ### Times wrt to DS horizontal average
        ToFTWt_rel=self.TDS0-ToFTWtime
        
        ### Make histograms
        dtvqdc=f'dtvqdc_{fixed_ch}_corrected'
        if not dtvqdc in hists:
            subtitle='{Time-walk corrected t_{0}^{DS}-t^{tw corr}_{SiPM '+str(SiPM)+'} v QDC_{SiPM '+str(SiPM)+'} w/ run'+str(self.TWCorrectionRun)+'};QDC_{SiPM '+str(SiPM)+'} [a.u];t_{0}^{DS}-t^{tw corr}_{SiPM '+str(SiPM)+'} [ns]'
            title='#splitline{'+ReadableFixedCh+'}'+subtitle
            hists[dtvqdc]=ROOT.TH2F(dtvqdc, title, 200, 0, 200, 800, -20, 20)         
      
        SiPMtimeToFcorrected=f'tSiPMToFcorrected_{fixed_ch}_corrected'
        if not SiPMtimeToFcorrected in hists:
            subtitle='{Time-walk and signal time-of-flight corrected SiPM time w/ run'+str(self.TWCorrectionRun)+'};t^{tw corr}_{SiPM '+str(SiPM)+'} [ns];Counts'
            title='#splitline{'+ReadableFixedCh+'}'+subtitle
            # hists[SiPMtimeToFTWcorrected]=ROOT.TH1F(SiPMtimeToFTWcorrected,title, 400, -10, 10)    
            if self.timealignment=='old': hists[SiPMtimeToFcorrected]=ROOT.TH1F(SiPMtimeToFcorrected, title, 600, -10, 20)
            else: hists[SiPMtimeToFcorrected]=ROOT.TH1F(SiPMtimeToFcorrected, title, 400, -10, 10)

        ### Fill histograms 
        hists[dtvqdc].Fill(qdc, ToFTWt_rel)
        hists[SiPMtimeToFcorrected].Fill(ToFTWtime)

    def WriteOutHistograms(self):                    

        for h in self.hists:

            if len(h.split('_'))==4:
                if len(h.split('_'))==4: histkey,detID,SiPM,state=h.split('_')
                fixed_ch='_'.join((detID,SiPM))
                hist=self.M.h[h]
                outpath=f'{self.outpath}/splitfiles/run{self.runNr}/{fixed_ch}/'
                path_obj=Path(outpath)
                path_obj.mkdir(parents=True, exist_ok=True)
                
                outfile=f'timewalk_{fixed_ch}_{self.options.nStart}.root'
                if os.path.exists(outpath+outfile): f=ROOT.TFile.Open(outpath+outfile, 'update')
                else: f=ROOT.TFile.Open(outpath+outfile, 'create')
                f.WriteObject(hist, hist.GetName(), 'kOverwrite')
                if self.options.mode=='zeroth':
                    f.WriteObject(self.hists['TDS0'], 'TDS0', 'kOverwrite')
                f.Close()
        print(f'{len(self.M.h)} histograms saved to {self.outpath}splitfiles/run{self.runNr}/fixed_ch')

    def LoadSystemObservables(self):
        systemobservablesfilename=f'{self.afswork}SystemObservables/run{self.TWCorrectionRun}/systemparameters.json'
        with open(systemobservablesfilename, 'r') as f:
            self.systemobservables=json.load(f)

    def yresidual3(self, detID):
        self.MuFilter.GetPosition(detID,self.A,self.B)
        doca=self.muAna.Getyresidual(detID)
        s,p,b=self.muAna.parseDetID(detID)
        if s==3: return True
        key=10*s+p
        
        dy_min, dy_max = self.dycut(self.cutdists[f'dy_{key}'])
        if doca>dy_max or doca<dy_min: return False
        else: return True

    def nSiPMscut(self, hit, nLeft, nRight):
        s,p,b=self.muAna.parseDetID(hit.GetDetectorID())
        if s==1: 
            if nLeft<6 or nRight<6: return False
        elif s==2:
            if nLeft<4 or nRight<4: return False
        elif s==3: 
            pass 
        return True   

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
            if td < mean-3*stddev or td > mean+3*stddev: return False
            else: return True

    def dycut(self, hist, nsig=1):    
        dymin=hist.GetMean()-nsig*hist.GetStdDev()
        dymax=hist.GetMean()+nsig*hist.GetStdDev()
        return dymin, dymax

    def slopecut(self, slopecut=0.1):
        slopeX, slopeY = self.mom.x()/self.mom.z(), self.mom.y()/self.mom.z()
        if abs(slopeX)>slopecut or abs(slopeY)>slopecut: return 0
        else: return 1

    def TDS0cut(self):
        if self.TDS0==-999. or self.TDS0==-998. or self.TDS0==-6237.5: return 0
        TDS0ns=self.TDS0*self.TDC2ns
        if TDS0ns<13.75 or TDS0ns>15.45: return 0
        return TDS0ns 
    
    def GetDistanceToSiPM(self,Ex):
        return ROOT.TMath.Sqrt((self.A.x()-Ex.x())**2+(self.A.y()-Ex.y())**2+(self.A.z()-Ex.z())**2)
