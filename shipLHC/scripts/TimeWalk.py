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

class TimeWalk(ROOT.FairTask):

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
        self.timingalignment=self.muAna.GetTimeAlignmentType(runNr=self.runNr)

        self.freq=160.316E6
        self.TDC2ns=1E9/self.freq

        self.largeSiPMmap={0:0 ,1:1 ,3:2 ,4:3 ,6:4 ,7:5}
        self.verticalBarDict={0:1, 1:3, 2:5, 3:6}
        self.xref=42. # To be set as USbarlength/2.
        
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
        
        if options.mode=='TW':
            if not options.CorrectionType: self.CorrectionType=4

            self.correctionparams=lambda ps : [y for x,y in enumerate(ps) if x%2==0]
            self.CorrectionType=options.CorrectionType
            if options.CorrectionType==1: self.correctionfunction = lambda ps, qdc : 1/sum( [ ps[i]*qdc**i for i in range(len(ps)) ] ) 
            elif options.CorrectionType==5: self.correctionfunction = lambda ps, qdc : ps[0]*ROOT.Log(ps[1]*qdc)+ps[2]
            elif options.CorrectionType==4: self.correctionfunction = lambda ps, qdc : ps[3]*(qdc-ps[0])/( ps[1] + ps[2]*(qdc-ps[0])*(qdc-ps[0]) ) 

            if options.TWCorrectionRun==None:
                if self.timingalignment=='old': self.TWCorrectionRun=str(5097).zfill(6)
                elif self.timingalignment=='new': self.TWCorrectionRun=str(5408).zfill(6)
            else:
                if self.muAna.GetTimeAlignmentType(runNr=options.TWCorrectionRun) != self.timingalignment:
                    if self.timingalignment=='old': self.TWCorrectionRun=str(5097).zfill(6)
                    elif self.timingalignment=='new': self.TWCorrectionRun=str(5408).zfill(6)
                else: self.TWCorrectionRun=str(options.TWCorrectionRun).zfill(6)

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

        if self.options.debug:
            self.FillChannelRateHists()

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
        # slopes=[i:(momvectors[i].x()/momvectors[i].z(), momvectors[i].y()/momvectors[i].z()) for i in momvectors]
        # if not TimeWalk.slopecut(momvectors[3]): return
        if not TimeWalk.slopecut(self.mom): return
        
        #### Testing clusters
        if self.options.debug:
            print(self.M.EventNumber)
        
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

            pred=self.GetDistanceToSiPM(A,Ex)

            channels_t=hit.GetAllTimes()
            channels_qdc=hit.GetAllSignals()

            for channel in channels_t:
                SiPM,clock=channel
                qdc=self.muAna.GetChannelVal(SiPM, channels_qdc)
                if qdc==-999.: continue
                fixed_ch=f'{detID}_{SiPM}'
                if self.options.mode=='zeroth': self.zeroth(fixed_ch, pred, clock, qdc, TDS0)
                elif self.options.mode=='ToF': self.ToF(fixed_ch, pred, clock, qdc, TDS0)
                elif self.options.mode=='TW': self.TW(fixed_ch, pred, clock, qdc, TDS0, meantimecorrection=False)
                elif self.options.mode=='globalcorrection':self.TW(fixed_ch, pred, clock, qdc, TDS0, meantimecorrection=True)

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

    def zeroth(self,fixed_ch,pred,clock,qdc,TDS0):
        hists=self.hists
        
        detID, SiPM=( int(fixed_ch.split('_')[i]) for i in range(2) ) 
        s,p,b=self.muAna.parseDetID(detID)

        correctedtime=clock*self.TDC2ns
        t_rel=TDS0-correctedtime
        if not self.muAna.DSVcheck(detID):  dtvpred=f'dtvxpred_{fixed_ch}_{self.state}'
        else: dtvpred=f'dtvypred_{fixed_ch}_{self.state}'
        
        ReadableDetID=self.muAna.MakeHumanReadableDetID(fixed_ch)
    
        if not dtvpred in hists:
            coord='x' if not self.muAna.DSVcheck(detID) else 'y'
            title='Uncorrected t_{0}^{DS}-t_{SiPM} v '+coord+'-position'
            splittitle='#splitline{'+ReadableDetID+'}{'+title+'}'
            axestitles=coord+'_{predicted} [cm];t_{0}^{DS}-t_{SiPM} [ns]'
            fulltitle=splittitle+';'+axestitles
            hists[dtvpred]=ROOT.TH2F(dtvpred,fulltitle,110,-10,100, 800, -20, 20.)

        SiPMtime=f'tSiPM_{fixed_ch}_{self.state}'
        if not SiPMtime in hists:
            title='Uncorrected t_{SiPM}'
            splittitle='#splitline{'+ReadableDetID+'}{'+title+'}'
            axestitles='t_{SiPM} [ns];Counts'
            fulltitle=splittitle+';'+axestitles            
            hists[SiPMtime]=ROOT.TH1F(SiPMtime,fulltitle, 400, -10, 10)
        
        attlen=f'attlen_{fixed_ch}_{self.state}'
        if not attlen in hists:
            coord='x' if not self.muAna.DSVcheck(detID) else 'y'
            title='Predicted position against QDC_{SiPM}'
            splittitle='#splitline{'+ReadableDetID+'}{'+title+'}'
            axestitles=coord+'_{predicted} [cm];QDC_{SiPM} [a.u]'
            fulltitle=splittitle+';'+axestitles 
            hists[attlen]=ROOT.TH2F(attlen,fulltitle, 110, -10, 110, 200, 0., 200)

        if self.options.debug and s==3 and self.muAna.DSVcheck(detID):
            verticalplanehitrate=f'verticalplanehitrate_plane{p}'
            if not verticalplanehitrate in hists:
                title=f'Vertical plane hit rate for plane {p+1};Bar number;Counts'
                hists[verticalplanehitrate]=ROOT.TH1F(verticalplanehitrate,title, 120, 0, 119)
            hists[verticalplanehitrate].Fill(b)
        
        self.hists[dtvpred].Fill(pred,t_rel)
        self.hists[SiPMtime].Fill(correctedtime)
        self.hists[attlen].Fill(pred, qdc)

    def ToF(self, fixed_ch, pred, clock, qdc, TDS0):        
        hists=self.hists
        ReadableDetID=self.muAna.MakeHumanReadableDetID(fixed_ch)

        cdata=self.muAna.Getcscint(self.runNr, fixed_ch, self.state)
        if cdata==-999.: return 0
        SiPM=int(fixed_ch.split('_')[-1])
        correctedtime=self.muAna.correct_ToF(SiPM, clock, pred, cdata, self.xref)[1]
        dtvqdc=f'dtvqdc_{fixed_ch}_{self.state}'
        if not dtvqdc in hists:
            subtitle='Uncorrected t_{0}^{DS}-t_{SiPM} v QDC_{SiPM};QDC_{SiPM} [a.u];t_{0}^{DS}-t_{SiPM} [ns]'
            title='#splitline{'+ReadableDetID+'}{'+subtitle+'}'
            hists[dtvqdc]=ROOT.TH2F(dtvqdc, title, 200, 0, 200, 800, -20, 20)
        t_rel=TDS0-correctedtime
        self.hists[dtvqdc].Fill(qdc,t_rel)

    def TW(self, fixed_ch, pred, clock, qdc, TDS0, meantimecorrection=False):
        hists=self.hists
        ReadableDetID=self.muAna.MakeHumanReadableDetID(fixed_ch)
        
        detID=int(fixed_ch.split('_')[0])
        SiPM=int(fixed_ch.split('_')[-1])
        time=clock*self.TDC2ns
        
        tmp = self.muAna.GetPolyParams(self.TWCorrectionRun, fixed_ch, 'uncorrected', self.CorrectionType)
        if tmp==-999.:
            return
        paramsAndErrors, correctionlimit, correctedmeantime=tmp
        if not correctedmeantime:
            # print(f'No corrected mean time for {fixed_ch}')
            meantimecorrection=False
        if meantimecorrection: correctedmeantime=float(correctedmeantime)

        chi2pNDF = self.muAna.Getchi2pNDF(self.TWCorrectionRun, fixed_ch, self.CorrectionType, 'uncorrected')
        cdata = self.muAna.Getcscint(self.TWCorrectionRun, fixed_ch, 'uncorrected')
        
        if paramsAndErrors==-999. or cdata==-999.: return 0
        polyparams = self.correctionparams(paramsAndErrors)
       
        ToFcorrectedtime=self.muAna.correct_ToF(SiPM, clock, pred, cdata, self.xref)[1]
        polycorrection = self.correctionfunction(polyparams, qdc)

        ### TW corrected time then ToF & TW corrected time        
        TWcorrectedtime=time+polycorrection
        ToFTWcorrectedtime=ToFcorrectedtime+polycorrection

        ### Times wrt to DS horizontal average
        ToFTWt_rel=TDS0-ToFTWcorrectedtime
        TWt_rel=TDS0-TWcorrectedtime

        ### Make histograms
        if not self.muAna.DSVcheck(detID):  dtvxpred=f'dtvxpred_{fixed_ch}_{self.state}'
        else: dtvxpred=f'dtvypred_{fixed_ch}_{self.state}'
        if not dtvxpred in hists:
            coord='x' if not self.muAna.DSVcheck(detID) else 'y'
            subtitle='Corrected t_{0}^{DS}-t_{SiPM} v '+coord+'-position w/ run'+str(self.TWCorrectionRun)+';'+coord+'_{predicted} [cm];t_{0}^{DS}-t_{SiPM} [ns]'
            title='#splitline{'+ReadableDetID+'}{'+subtitle+'}'
            hists[dtvxpred]=ROOT.TH2F(dtvxpred,title,110,-10,100, 800, -20, 20.)
        
        # dtvqdc=f'dtvqdc_{fixed_ch}_{self.state}'
        dtvqdc=f'dtvqdc_{fixed_ch}_corrected'
        if not dtvqdc in hists:
            subtitle='Corrected t_{0}^{DS}-t_{SiPM} v QDC_{SiPM} w/ run'+str(self.TWCorrectionRun)+';QDC_{SiPM} [a.u];t_{0}^{DS}-t_{SiPM} [ns]'
            title='#splitline{'+ReadableDetID+'}{'+subtitle+'}'
            hists[dtvqdc]=ROOT.TH2F(dtvqdc, title, 200, 0, 200, 800, -20, 20)                               

        SiPMtime=f'tSiPM_{fixed_ch}_corrected'
        if not SiPMtime in hists:
            subtitle=f'Corrected SiPM time w/ run'+str(self.TWCorrectionRun)+';SiPM time [ns];Counts'
            title='#splitline{'+ReadableDetID+'}{'+subtitle+'}'
            hists[SiPMtime]=ROOT.TH1F(SiPMtime,title, 400, -10, 10)
        
        ### Fill histograms 
        hists[dtvxpred].Fill(pred, TWt_rel)
        hists[dtvqdc].Fill(qdc, ToFTWt_rel)
        hists[SiPMtime].Fill(TWt_rel)
        
        if meantimecorrection:

            tds0correctedtwcorrectedtime=TWcorrectedtime-correctedmeantime

            dtvqdc_meantimecorrected=f'dtvqdc_{fixed_ch}_corrected_tds0corrected'
            if not dtvqdc_meantimecorrected in hists:
                subtitle='TW and t_{0}^{DS}-mean corrected t_{0}^{DS}-t_{SiPM} v QDC_{SiPM} w/ run'+str(self.TWCorrectionRun)+';QDC_{SiPM} [a.u];t_{0}^{DS}-t_{SiPM} [ns]'
                title='#splitline{'+ReadableDetID+'}{'+subtitle+'}'
                hists[dtvqdc_meantimecorrected]=ROOT.TH2F(dtvqdc_meantimecorrected, title, 200, 0, 200, 800, -20, 20)   

            SiPMtime_meantimecorrected=f'tSiPM_{fixed_ch}_corrected_meantimecorrected'
            if not SiPMtime_meantimecorrected in hists:
                subtitle='TW and t_{0}^{DS}-mean corrected SiPM time w/ run'+str(self.TWCorrectionRun)+';SiPM time [ns];Counts'
                title='#splitline{'+ReadableDetID+'}{'+subtitle+'}'
                hists[SiPMtime_meantimecorrected]=ROOT.TH1F(SiPMtime_meantimecorrected,title, 500, 0, 25)
            
            hists[dtvqdc_meantimecorrected].Fill(qdc, tds0correctedtwcorrectedtime)
            hists[SiPMtime_meantimecorrected].Fill(tds0correctedtwcorrectedtime)

    def AllCorrection(self, fixed_ch, pred, clock, qdc, TDS0, hit):
        pass

        # if self.options.debug==1: hists[digimethod].Fill(qdc, DigiTWToFt_rel)

    def WriteOutHistograms(self):
        testing_outfile=f'testing_{self.options.runNumber}_{self.options.nStart}_{self.options.nEvents}.root'
        testing_f=ROOT.TFile.Open(testing_outfile, 'update')                    

        for h in self.M.h:
            if self.options.debug:
                if h.find('channelhitrate')!=-1:
                    hist=self.M.h[h]
                    testing_f.WriteObject(hist, f'TimeWalk-{h}')
            if len(h.split('_'))==4 or any( [h.find(i)!=-1 for i in ('tds0corrected', 'meantimecorrected')] ):
                if len(h.split('_'))==4: histkey,detID,SiPM,state=h.split('_')
                elif any( [h.find(i)!=-1 for i in ('tds0corrected', 'meantimecorrected')] ):
                    histkey, detID, SiPM, state, extrakey=h.split('_')
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
        
        if len(US1hits)==0: 
            return -999
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
        if self.timingalignment=='old':
            if td<0: return False 
            else: return True
        else: 
            if not 'timingdiscriminant' in self.cutdists:
                print(f'No timing discriminant histogram')
                return 
            hist=self.cutdists['timingdiscriminant']
            mean, stddev=hist.GetMean(), hist.GetStdDev()
            # print(f'td: {td}, mean +/- 2*std. dev.: {mean-2*stddev}, {mean+2*stddev}')
            if td < mean-3*stddev or td > mean+3*stddev: return False
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
        return ROOT.TMath.Sqrt((A.x()-Ex.x())**2+(A.y()-Ex.y())**2+(A.z()-Ex.z())**2)

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
                      