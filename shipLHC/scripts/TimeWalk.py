import ROOT, os, csv, json
from argparse import ArgumentParser
from AnalysisFunctions import Analysis
from pathlib import Path
import numpy as np

deciles=[i/10 for i in range(11)]
SiPM_distance_cut=5

ROOT.gStyle.SetTitleXOffset(0.7)
ROOT.gStyle.SetTitleYOffset(0.6)
for i in ('x', 'y', 'z','t'): ROOT.gStyle.SetTitleSize(0.04, i)

class TimeWalk(ROOT.FairTask):

    def __init__(self, options, monitor):

        self.M=monitor 
        self.simulation = options.simulation
        self.muAna = Analysis(options)
        self.options=options

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
        self.A, self.B = ROOT.TVector3(), ROOT.TVector3()  

        freq=160.316E6
        self.TDC2ns=1E9/freq

        self.largeSiPMmap={0:0 ,1:1 ,3:2 ,4:3 ,6:4 ,7:5}
        self.verticalBarDict={0:1, 1:3, 2:5, 3:6}        
     
        self.hists=self.M.h

        run=ROOT.FairRunAna.Instance()

        ioman=ROOT.FairRootManager.Instance()
        self.OT=ioman.GetSink().GetOutTree()

        """
        NEW FEATURE!
        
        In order to switch between analysis with different reference systems, now I can submit runs with mode {mode}{referencesystem} and parse the options.mode here and look for an integer, if an integer is found that that will be passed as 
        options.referencesystem. By default options.referencesystem is 3 (DS)
        """
        if options.mode.find('1')>0: 
            self.referencesystem=1
            self.mode=options.mode[:-1]
        elif options.mode.find('3')>0: 
            self.referencesystem=3
            self.mode=options.mode[:-1]
        else:
            self.referencesystem=options.referencesystem
            self.mode=options.mode
        self.refsysname='DS' if self.referencesystem==3 else 'SF'   

        if options.path.find('commissioning/TI18')>0:
            self.outpath=options.afswork+'-commissioning/'
            self.path='TI18'
        elif options.path.find('physics/202')>0: 
            self.outpath=options.afswork+'-physics2022/'
            self.path='TI18'
        elif options.path.find('testbeam_June2023_H8')>0: 
            self.outpath=options.afswork+'-H8/'
            self.path='H8'
        elif options.path=='./':
            self.outpath=options.afswork+'-physics2022/'
            self.path='TI18'

        #### Time walk correction run is independent of time alignment settings. 
        self.TWCorrectionRun=str(5408).zfill(6)            
        
        if self.simulation:
            self.timealignment='sim'
            self.GetLumis() # loads the simulation luminosities dictionary into tw.simlumis
            self.simMode = options.simMode
            self.cutdists=self.muAna.GetCutDistributions("005408", ('dy', 'timingdiscriminant'))

        elif not self.simulation:
            statedict={'zeroth':'uncorrected', 'tof':'uncorrected',
                        'tw':'corrected', 'res':'corrected', 
                        'selectioncriteria':'corrected',
                        'systemalignment':'corrected',
                        'reconstructmuonposition':'corrected',
                        'showerprofiles':'corrected', 
                        'numusignalevents':'corrected',
                        'tds0-studies':'uncorrected',
                        'extendedreconstruction':'corrected', 'struckquark':'corrected'}
            
            self.state=statedict[self.mode]
            
            self.timealignment=self.muAna.GetTimeAlignmentType(runNr=self.runNr)
            
            if options.debug: self.trackevents=[]

            ### If no time-walk correction run is provided. Set the default correction run depending on time alignment of the data set
            if options.AlignmentRun==None or self.muAna.GetTimeAlignmentType(runNr=options.AlignmentRun) != self.timealignment:
                if self.timealignment=='old': self.AlignmentRun=str(5097).zfill(6)
                elif self.timealignment=='new': self.AlignmentRun=str(5408).zfill(6)
                elif self.timealignment=='new+LHCsynch': self.AlignmentRun=str(5999).zfill(6)
            self.muAna.AlignmentRun=self.AlignmentRun

            self.cutdists=self.muAna.GetCutDistributions(self.TWCorrectionRun, ('dy', 'timingdiscriminant'))

            self.muAna.MakeAlignmentParameterDict()
            self.muAna.Makecscintdict(self.TWCorrectionRun, state=self.state)
            self.muAna.MakeTWCorrectionDict()

        if self.mode == 'systemalignment':
            from systemalignment import SystemAlignment
            self.sa = SystemAlignment(options, self)

        if self.mode in ['reconstructmuonposition', 'multimuonreconstruction']:
            from systemalignment import SystemAlignment
            self.sa = SystemAlignment(options, self)  

        elif self.mode == 'showerprofiles':
            from showerprofiles import ShowerProfiles
            if options.numuStudy: self.numuStudy=True
            else: self.numuStudy=False
            self.sp = ShowerProfiles(options, self)
            self.barycentredata={}

        elif self.mode == 'selectioncriteria':
            from selectioncriteria import MuonSelectionCriteria as SelectionCriteria
            self.sc = SelectionCriteria(options, self)   
            # self.sc.cuteffect['no cuts'] = self.M.eventTree.GetEntries()

        elif self.mode.find('extendedreconstruction')>-1:
            from extendedmuonreconstruction import ExtendedMuonReconstruction as ExtendedMuonReconstruction
            self.emr = ExtendedMuonReconstruction(options, self)
        elif self.mode=='struckquark':
            from extendedmuonreconstruction import QuarkVectorExtrapolation as QuarkVectorExtrapolation
            self.qve = QuarkVectorExtrapolation(options, self)

        with open(f'/afs/cern.ch/user/a/aconsnd/Timing/TWhistogramformatting.json', 'r') as x:
            self.histformatting=json.load(x)  

        self.notInDS=0

        self.trackRequired=True if self.mode in('zeroth', 'tof', 'tw', 'res', 'systemalignment', 'reconstructmuonposition', 'selectioncriteria') else False

    def GetEntries(self):
        return self.eventTree.GetEntries()

    def InvestigateSignalEvents(self):
        self.numutiming.InvestigateSignalEvents()

    def GetEventScaledWeight(self):
        w = self.M.eventTree.GetWeight()
        luminosity = self.simlumis[self.simMode] 
        return round(1/(w*luminosity), 4)

    def GetLumis(self):
        simlumi_filename = f'/afs/cern.ch/user/a/aconsnd/Timing/simulationluminosities.json'
        with open(simlumi_filename, 'r') as x:
            self.simlumis = json.load(x)

    def ExecuteEvent(self, event):

        if not hasattr(self.muAna, 'task'): self.muAna.SetTask(self)

        tracks={1:[], 3:[]}
        inScifi, inDS = False, False

        hits = event.Digi_MuFilterHits
        scifi_hits = event.Digi_ScifiHits

        if self.simulation: 
            self.scaled_event_weight = self.GetEventScaledWeight()
            histname = f'{self.simMode}_eventweight'
            if not histname in self.hists:
                title = f'Event weight for: {self.simMode};Event weight;Counts'
                self.hists[histname] = ROOT.TH1D(histname, title, 200, -1E3, 1E3)
            self.hists[histname].Fill(self.M.eventTree.GetWeight())

        for i,track in enumerate(self.M.Reco_MuonTracks):
            if any([not track.getFitStatus().isFitConverged(), track.getFitStatus().getNdf()==0]): 
                continue
            if track.GetUniqueID()==1:
                inScifi=True
                tracks[1].append(self.M.Reco_MuonTracks[i])
            if track.GetUniqueID()==3: 
                inDS=True
                tracks[3].append(self.M.Reco_MuonTracks[i])

        ### If there are more than 1 DS track, take the track with the lowest chi2/Ndf
        if len(tracks[self.referencesystem])==1: self.track = tracks[self.referencesystem][0]
        elif len(tracks[self.referencesystem])>1: 
            tmp={i:i.getFitStatus().getChi2()/i.getFitStatus().getNdf() for i in tracks[self.referencesystem]}
            self.track=tmp[ min(tmp) ]

        self.hasTrack=True
        if (self.referencesystem==3 and not inDS) or (self.referencesystem==1 and not inScifi): 
            self.hasTrack=False 
        
        # Only throw event out if we need a track but don't have one
        if self.trackRequired and not self.hasTrack: return 

        elif self.hasTrack:       

            # cut flow hist
            if self.mode=='selectioncriteria': self.sc.cuteffect['has track']+=1

            fitStatus=self.track.getFitStatus()
            fstate=self.track.getFittedState()

            self.pos=fstate.getPos()
            self.mom=fstate.getMom()
            self.trackchi2NDF=fitStatus.getChi2()/fitStatus.getNdf()+1E-10
            self.chi2pndf_hist(fitStatus)

            # cut flow hist
            self.passredchi2cut = False 
            if self.trackchi2NDF < 5: self.passredchi2cut = True
            if self.mode=='selectioncriteria' and self.passredchi2cut and self.hasTrack: self.sc.cuteffect['#chi^{2}_{#nu}']+=1

            # Get DS event t0
            if self.referencesystem==3 and self.hasTrack and not self.simulation:
                x=self.muAna.GetDSHaverage(hits) # Now returns in ns
                if x==-999.: return
                self.reft, firedDSHbars=x
                if not 'reft-DS' in self.hists:
                    self.hists['reft-DS']=ROOT.TH1F('reft','Average time of DS horizontal bars;DS horizontal average time [ns];Counts', 100, 0, 25)
                self.hists['reft-DS'].Fill(self.reft) 

            elif self.referencesystem==1 and self.hasTrack and not self.simulation:
                self.reft = self.muAna.GetScifiAverageTime(self.Scifi, scifi_hits)
                if not self.reft: return
                if not 'reft-SF' in self.hists:
                    self.hists['reft-SF']=ROOT.TH1F('reft','Average time of SiPMs in Scifi;Time [ns];Counts', 200, -25, 25)
                self.hists['reft-SF'].Fill(self.reft)

            # self.passtds0cut = self.tds0cut()
            # if self.mode=='selectioncriteria' and self.passtds0cut and self.passredchi2cut and self.hasTrack: 
            #     self.sc.cuteffect['t^{ds}_{0}']+=1

            self.Ex = {}
            for plane in range(5):
                zEx=self.zPos['MuFilter'][20+plane]
                lam=(zEx-self.pos.z())/self.mom.z()
                self.Ex[plane]=ROOT.TVector3(self.pos.x()+lam*self.mom.x(), self.pos.y()+lam*self.mom.y(), self.pos.z()+lam*self.mom.z())
            
            self.passslopecut = self.slopecut()
            if self.mode=='selectioncriteria' and self.passslopecut and self.passredchi2cut and self.hasTrack:
                self.sc.cuteffect['slope']+=1  

        # Get DS event t0
        if self.referencesystem==3 and self.hasTrack and not self.simulation:

            t0vxpred=f'tds0vxpred'
            if not t0vxpred in self.hists:
                title='t_{0}^{DS} v x-position expected at HCAL plane 5;x-position at plane 5 [cm];t_{DS}^{0}'
                self.hists[t0vxpred]=ROOT.TH2F(t0vxpred, title, 110, 10, -100, 200, 0, 50) 
            self.hists[t0vxpred].Fill(self.Ex[4].x(), self.reft)

            self.td = self.muAna.GetTimingDiscriminant(hits) # Require that US1 TDC average is less than the DSH TDC average to ensure forward travelling track
            if self.TimingDiscriminantCut(): self.passtdcut=True 
            else: self.passtdcut=False

            # if self.mode=='selectioncriteria' and self.passtdcut and self.passslopecut and self.passredchi2cut and self.hasTrack:
            #     self.sc.cuteffect['timingdiscriminant']+=1
        
        if self.mode == 'selectioncriteria':
            self.sc.Filltds0vcuts()
            self.sc.FillHists(hits, scifi_hits)
            return
        
        elif self.mode=='showerprofiles':
            # self.sp.ExtractScifiData(scifi_hits)
            self.sp.ShowerDirection(hits, scifi_hits)
            self.sp.EvaluateShowerSlopes() # Quite slow, comment is not needed
            return
                
        elif self.mode.find('extendedreconstruction')>-1:
            self.emr.ExtendReconstruction(hits, scifi_hits)
            return
        
        elif self.mode == 'struckquark': 
            self.qve.StruckQuarkExtrapolation(hits)
            return
        
        elif self.mode == 'struckquark': 
            self.qve.StruckQuarkExtrapolation(hits)
            return
        
        elif self.mode == 'reconstructmuonposition':
            print(f'running rmp code')
            self.sa.ReconstructMuonPosition(hits)
            return        

        ### Slope cut
        if self.slopecut(): self.passslopecut=True 
        else: self.passslopecut=False

        ### Slope cut
        if self.slopecut(): self.passslopecut=True 
        else: self.passslopecut=False

        for hit in event.Digi_MuFilterHits:
            nLeft, nRight=self.muAna.GetnFiredSiPMs(hit)
            if not hit.isValid(): continue

            detID=hit.GetDetectorID()
            s,p,b=self.muAna.parseDetID(detID)            
            if s!=2: continue

            # Only investigate veto hits if there is a Scifi track!
            if not inScifi and s==1: continue

            self.MuFilter.GetPosition(detID,self.A,self.B)

            if self.yresidual3(detID): self.trackrelated = True
            else: self.trackrelated=False

            if self.options.debug: self.trackevents.append(self.M.EventNumber)

            channels_t=hit.GetAllTimes()
            channels_qdc=hit.GetAllSignals()

            # Only investigate track related hits
            if not self.trackrelated:continue             
            
            if not all([self.passslopecut, self.passredchi2cut]): return 

            if self.mode=='systemalignment' and s==2:
                self.sa.FillSiPMHists(hit)
                if self.options.XT: self.sa.XTHists(hit)
                self.sa.FillBarHists(hit)
                # self.sa.ScifiCorrectedTimes(hit)
                continue

            # Hist just used for discussion in the thesis
            if self.mode=='zeroth': self.MakeExplantoryHist(event.Digi_MuFilterHits)

            for channel in channels_t:
                SiPM,clock=channel
                qdc=self.muAna.GetChannelVal(SiPM, channels_qdc)
                if qdc==-999.: continue
                fixed_ch=f'{detID}_{SiPM}'
                if self.mode=='zeroth': self.zeroth(fixed_ch, clock, qdc)
                elif self.mode=='tof': self.tof(fixed_ch, clock, qdc)
                elif self.mode=='tw': self.tw(fixed_ch, clock, qdc, meantimecorrection=False)
                elif self.mode=='res': self.res(fixed_ch, clock, qdc)

    def chi2pndf_hist(self, fitStatus):
        histname = 'reducedchi2'
        if not histname in self.hists:
            if self.referencesystem==3: tmp='muon system'
            else: tmp='scifi'
            title = '#chi^{2} per degree of freedom for tracks in '+tmp+';#chi^{2}_{#nu} [unitless];Counts'
            self.hists[histname] = ROOT.TH1F(histname, title, 100, 0, 50)
        self.hists[histname].Fill(self.trackchi2NDF)

        histname = 'chi2'
        if not histname in self.hists:
            if self.referencesystem==3: tmp='muon system'
            else: tmp='scifi'
            title = '#chi^{2} for tracks in '+tmp+';#chi^{2} [unitless];Counts'
            self.hists[histname] = ROOT.TH1F(histname, title, 100, 0, 100)
        self.hists[histname].Fill(fitStatus.getChi2())

        histname = 'dof'
        if not histname in self.hists:
            if self.referencesystem==3: tmp='muon system'
            else: tmp='scifi'
            title = 'degrees of freedom for tracks in '+tmp+';Degrees of freedom [unitless];Counts'
            self.hists[histname] = ROOT.TH1F(histname, title, 5, 0, 5)
        self.hists[histname].Fill(fitStatus.getNdf())

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

    def fiducialvolumecut(self, fixed_ch):
        detID, ch = fixed_ch.split('_')

        if int(ch)>=8: side='right'
        else: side='left'

        s,p,b = self.muAna.parseDetID(detID)
        self.MuFilter.GetPosition(int(detID), self.A, self.B)

        # Do not accept hits where the muon extrapolates to within 2 cm of the bar end
        if side == 'left' and (self.A.x() - self.Ex[p].x()) < SiPM_distance_cut: 
            print(f'Bar end: {self.A.x()}, Ex: {self.Ex[p].x()}')
            return False 
        elif side == 'right' and (self.Ex[p].x() - self.B.x()) < SiPM_distance_cut:
            print(f'Bar end: {self.B.x()}, Ex: {self.Ex[p].x()}')
            return False 
        else: return True

    def zeroth(self,fixed_ch,clock,qdc):

        # if not all([self.passslopecut, self.passredchi2cut]): return 

        hists=self.hists
        
        detID, SiPM = ( int(fixed_ch.split('_')[i]) for i in range(2) )
        s,p,b=self.muAna.parseDetID(detID)

        correctedtime=clock*self.TDC2ns
        t_rel=self.reft-correctedtime
        if not self.muAna.DSVcheck(detID):  dtvpred=f'dtvxpred_{fixed_ch}_{self.state}'
        else: dtvpred=f'dtvypred_{fixed_ch}_{self.state}'
        
        ReadableFixedCh=self.muAna.MakeHumanReadableFixedCh(fixed_ch)
    
        if not dtvpred in hists:
            coord='x' if not self.muAna.DSVcheck(detID) else 'y'
            title='{No time-walk correction t_{0}^{DS}-t^{uncorr}_{SiPM} v '+coord+'-position}'
            splittitle='#splitline{'+ReadableFixedCh+'}'+title
            axestitles=coord+'_{predicted} [cm];t_{0}^{DS}-t^{uncorr}_{SiPM} [ns]'
            fulltitle=splittitle+';'+axestitles
            histformat = self.histformatting["dtvxpred"][self.state]
            hists[dtvpred]=ROOT.TH2F(dtvpred,fulltitle,*histformat[0], *histformat[1])

        attlen=f'attlen_{fixed_ch}_{self.state}'
        if not attlen in hists:
            coord='x' if not self.muAna.DSVcheck(detID) else 'y'
            title='{Predicted position against QDC_{SiPM}}'
            splittitle='#splitline{'+ReadableFixedCh+'}'+title
            axestitles=coord+'_{predicted} [cm];QDC_{SiPM} [a.u]'
            fulltitle=splittitle+';'+axestitles 
            hists[attlen]=ROOT.TH2F(attlen,fulltitle, 110, 10, -100, 200, 0., 200)

        self.hists[dtvpred].Fill(self.Ex[p].x(),t_rel)
        self.hists[attlen].Fill(self.Ex[p].x(), qdc)

    def MakeExplantoryHist(self, hits):
        # Make a QDC distribution for hits in each third of the bar        
        # v1, v2 = ROOT.TVector3(), ROOT.TVector3()

        for hit in hits: 
            detID=hit.GetDetectorID()
            # self.MuFilter.GetPosition(detID, v1, v2)
            s,p,b=self.muAna.parseDetID(detID)
            if not s==2: continue 
            
            for x in hit.GetAllSignals():
                SiPM, qdc = x

                if SiPM<8: distance_to_SiPM = self.A.x() - self.Ex[p].x()
                else: distance_to_SiPM = self.Ex[p].x() - self.B.x()
                
                bar_L = (self.B.x() - self.A.x())

                relative_bar_pos = (self.Ex[p].x() - self.B.x()) / bar_L # fraction of bar length
                if relative_bar_pos > 0 and relative_bar_pos <= 1/3: bar_third=1
                elif relative_bar_pos > 1/3 and relative_bar_pos <= 2/3: bar_third=2
                elif relative_bar_pos > 2/3 and relative_bar_pos <= 1: bar_third=3 
                else: 
                    # print(f'detID: {detID}, x_pred: {self.Ex[p].x()}')
                    continue 

                histname = f'qdcthird{bar_third}_{detID}_{SiPM}_{self.state}'
                if not histname in self.hists:
                    title=f'QDC distribution for hits {bar_third}/3 from the right side;QDC [a.u];Counts'
                    self.hists[histname] = ROOT.TH1F(histname, title, 200, 0, 200)
                self.hists[histname].Fill(qdc)
                
    def tof(self, fixed_ch, clock, qdc):

        # if not all([self.passslopecut, self.passredchi2cut]): return 

        hists=self.hists
        ReadableFixedCh=self.muAna.MakeHumanReadableFixedCh(fixed_ch)

        detID, SiPM=fixed_ch.split('_')
        s,p,b = self.muAna.parseDetID(int(detID))

        # if not fixed_ch in self.muAna.cscintvalues:return
        # cdata = self.muAna.cscintvalues[fixed_ch]
        # ToFcorrectedtime = self.muAna.correct_ToF(self.MuFilter, fixed_ch, clock, self.Ex[p].x())[1]
        # t_rel = self.reft - ToFcorrectedtime
        ToFcorrectedtime = self.muAna.MuFilterCorrectedTime(self.MuFilter, fixed_ch, qdc, clock, self.Ex[p].x())
        t_rel = self.reft - ToFcorrectedtime

        # Testing fiducial volume cut
        if self.fiducialvolumecut(fixed_ch): pass_fiducialvolumecut=True 
        else: pass_fiducialvolumecut = False

        dtvqdc = f'dtvqdc_{fixed_ch}_{self.state}'
        
        for cut in (None, f'-{SiPM_distance_cut}cmFiducialCut'):

            if cut!=None: histname = dtvqdc + cut
            else: histname = dtvqdc

            if not histname in hists:
                if cut==f'-{SiPM_distance_cut}cmFiducialCut': 
                    subtitle='{No tw correction t_{0}^{DS}-t^{uncorr}_{SiPM} v QDC_{SiPM}};QDC_{SiPM} [a.u];t_{0}^{DS}-t^{uncorr}_{SiPM} [ns]'
                else: subtitle='{'+str(SiPM_distance_cut)+' cm fiducial cut, no tw correction t_{0}^{DS}-t^{uncorr}_{SiPM} v QDC_{SiPM}};QDC_{SiPM} [a.u];t_{0}^{DS}-t^{uncorr}_{SiPM} [ns]'
                
                title='#splitline{'+ReadableFixedCh+'}'+subtitle
                histformat=self.histformatting["dtvqdc"][self.state]
                hists[histname]=ROOT.TH2F(histname, title, *histformat[0],*histformat[1])

        self.hists[dtvqdc].Fill(qdc,t_rel)
        if pass_fiducialvolumecut: 
            self.hists[dtvqdc+f'-{SiPM_distance_cut}cmFiducialCut'].Fill(qdc,t_rel)

    def tw(self, fixed_ch, clock, qdc, meantimecorrection=False):

        # if not all([self.passslopecut, self.passredchi2cut]): return 

        hists=self.hists
        ReadableFixedCh=self.muAna.MakeHumanReadableFixedCh(fixed_ch)
        
        detID, SiPM=int(fixed_ch.split('_')[0]), int(fixed_ch.split('_')[1])
        s,p,b = self.muAna.parseDetID(detID)

        # time=clock*self.TDC2ns
        
        # ### Check if fixed_ch has tw corr params and a cscint value
        # if any([fixed_ch not in self.muAna.twparameters, fixed_ch not in self.muAna.cscintvalues]): return 

        # twparams=self.muAna.twparameters[fixed_ch]
        # twcorrection = self.correctionfunction(twparams, qdc)

        # ### TW corrected time then ToF & TW corrected time
        # TWcorrectedtime=time+twcorrection

        # ### Times wrt to DS horizontal average
        # TWt_rel = self.reft-TWcorrectedtime
        TWcorrectedtime = self.muAna.MuFilterCorrectedTime(self.MuFilter, fixed_ch, qdc, clock, x=0)
        if not TWcorrectedtime: return
        TWt_rel = self.reft - TWcorrectedtime

        ### Make histograms
        dtvxpred=f'dtvxpred_{fixed_ch}_{self.state}'
        if not dtvxpred in hists and self.referencesystem==3:
            subtitle='{Time-walk corrected t_{0}^{DS}-t^{tw corr}_{SiPM} v x-position w/ run'+str(self.TWCorrectionRun)+'};x_{predicted} [cm];t_{0}^{DS}-t^{tw corr}_{SiPM} [ns]'
            title='#splitline{'+ReadableFixedCh+'}'+subtitle
            histformat = self.histformatting['dtvxpred'][self.state]        
            hists[dtvxpred]=ROOT.TH2F(dtvxpred,title,*histformat[0], *histformat[1])          
    
        hists[dtvxpred].Fill(self.Ex[p].x(), TWt_rel)
        
    def res(self, fixed_ch, clock, qdc):
        
        # if not all([self.passslopecut, self.passredchi2cut]): return 

        hists=self.hists
        ReadableFixedCh=self.muAna.MakeHumanReadableFixedCh(fixed_ch)
        detID, SiPM = int(fixed_ch.split('_')[0]), int(fixed_ch.split('_')[1])
        s,p,b=self.muAna.parseDetID(detID)

        # ## Check if fixed_ch has tw corr params and a cscint value
        # if any([fixed_ch not in self.muAna.twparameters, fixed_ch not in self.muAna.cscintvalues]): return 

        # cdata=self.muAna.cscintvalues[fixed_ch]
        # twparams=self.muAna.twparameters[fixed_ch]

        # ToFtime=self.muAna.correct_ToF(self.MuFilter, fixed_ch, clock, self.Ex[p].x())[1] # method returns (SiPM, ToF corrected time)
        # twcorrection = self.muAna.correctionfunction(twparams, qdc)

        # ### TW corrected time then ToF & TW corrected time
        # ToFTWtime = ToFtime + twcorrection

        # ### Times wrt to DS horizontal average
        # ToFTWt_rel = self.reft - ToFTWtime
        
        ToFTWtime = self.muAna.MuFilterCorrectedTime(self.MuFilter, fixed_ch, qdc, clock, x=self.Ex[p].x())
        if not ToFTWtime: return
        ToFTWt_rel = self.reft - ToFTWtime

        ### Make histograms
        if self.referencesystem==3: dtvqdc=f'dtvqdc_{fixed_ch}_corrected'
        elif self.referencesystem==1: dtvqdc=f'dtSFvqdc_{fixed_ch}_corrected'

        if not dtvqdc in hists:
            subtitle='{Time-walk corrected t_{0}^{'+self.refsysname+'}-t^{tw corr}_{SiPM '+str(SiPM)+'} v QDC_{SiPM '+str(SiPM)+'} w/ run'+str(self.TWCorrectionRun)+'};QDC_{SiPM '+str(SiPM)+'} [a.u];t_{0}^{'+self.refsysname+'}-t^{tw corr}_{SiPM '+str(SiPM)+'} [ns]'
            title='#splitline{'+ReadableFixedCh+'}'+subtitle
            histformat = self.histformatting['dtvqdc'][self.state]
            hists[dtvqdc]=ROOT.TH2F(dtvqdc, title, *histformat[0], *histformat[1])

        ### Fill histograms 
        hists[dtvqdc].Fill(qdc, ToFTWt_rel)

    def tds0_studies(self, hits):

        if self.passslopecut:
            name=f'tds0-slopecut'
            if not name in self.hists:
                title = 't_{0}^{ds} with slope cut;t_{0}^{ds} [ns];Counts'
                self.hists[name]=ROOT.TH1F(name, title, 200, 0, 50)
            self.hists[name].Fill(self.reft)
        
        if self.passtdcut:
            name=f'tds0-tdcut'
            if not name in self.hists:
                title = 't_{0}^{ds} with timing discriminant cut;t_{0}^{ds} [ns];Counts'
                self.hists[name]=ROOT.TH1F(name, title, 200, 0, 50)
            self.hists[name].Fill(self.reft)

        if self.passtdcut and self.slopecut:
            name=f'tds0-tdcut+slopecut'
            if not name in self.hists:
                title = 't_{0}^{ds} with timing discriminant cut and slope cut;t_{0}^{ds} [ns];Counts'
                self.hists[name]=ROOT.TH1F(name, title, 200, 0, 50)
            self.hists[name].Fill(self.reft)

        ### Apply slope cut to all these to minimise non-IP1 contribution
        if not self.passslopecut: return

        for mode in ('testing-tds0','testing-deltastations'):
            testing_value = self.muAna.GetDSHaverage(hits, mode)
            
            if mode.find('tds0')!=-1:

                if testing_value==-999.: continue
                tds0, firedDSHbars = testing_value

                if not mode in self.hists:
                    title = f'{mode.replace("-", " ")}; [ns];Counts'
                    bins = (200, 0, 50)
                    self.hists[mode] = ROOT.TH1F(mode, title, *bins)
                self.hists[mode].Fill(tds0)

            else: 
                if isinstance(testing_value, bool): continue
                bins = (100, -12.5, 12.5)
                for d in ('32', '21'):
                    name = f'{mode}{d}'
                    if not name in self.hists:
                        title = f'{mode.replace("-", " ")} {d}; [ns];Counts'
                        self.hists[name] = ROOT.TH1F(name, title, *bins)
                    self.hists[name].Fill(testing_value[f'delta{d}'])

    def Scifi_t0_studies(self, scifi_hits):
        name = f'scifi_averagetime'
        if not name in self.hists:
            title = f'{mode.replace("-", " ")} {d}; [ns];Counts'
            self.hists[name] = ROOT.TH1F(name, title, *bins)
        self.hists[name].Fill(testing_value[f'delta{d}'])  

    def GetSimEngine(self):
        simEngine = self.options.geoFile.split('.')[1].split('-')[0]
        return simEngine              

    def WriteOutHistograms(self):

        if self.mode in ('zeroth', 'tof', 'tw', 'res'):
            for h in self.hists:

                if not len(h.split('_'))==4: continue
                    
                if len(h.split('_'))==4: histkey,detID,SiPM,state=h.split('_')
                
                fixed_ch='_'.join((detID,SiPM))
                hist=self.hists[h]
                outpath=f'{self.outpath}/splitfiles/run{self.runNr}/{fixed_ch}/'
                path_obj=Path(outpath)
                path_obj.mkdir(parents=True, exist_ok=True)
                outfile=f'timewalk_{fixed_ch}_{self.options.nStart}.root'
                if os.path.exists(outpath+outfile): f=ROOT.TFile.Open(outpath+outfile, 'update')
                else: f=ROOT.TFile.Open(outpath+outfile, 'create')
                f.WriteObject(hist, hist.GetName(), 'kOverwrite')
                if self.mode=='res':
                    if self.referencesystem==1: name='reft-SF'
                    else: name='reft-DS'
                    f.WriteObject(self.hists[name], name, 'kOverwrite')
                    if self.refsysname=='DS': f.WriteObject(self.hists['tds0vxpred'], 'tds0vxpred', 'kOverwrite')
                    f.WriteObject(self.hists['reducedchi2'], 'reducedchi2', 'kOverwrite')  
                f.Close()

            print(f'{len(self.M.h)} histograms saved to {self.outpath}splitfiles/run{self.runNr}/fixed_ch/timewalk_fixed_ch_{self.options.nStart}.root')
        
        elif self.mode == 'systemalignment':
            self.sa.WriteOutHistograms()
            
        elif self.mode == 'showerprofiles':
            self.sp.WriteOutHistograms()      
            
        elif self.mode == 'selectioncriteria':
            self.sc.WriteOutHistograms()

        elif self.mode == 'reconstructmuonposition':
            self.sa.WriteOutReconstructionHistograms()   

        elif self.mode.find('extendedreconstruction')>-1:
            self.emr.WriteOutHistograms()

        elif self.mode.find('struckquark')>-1:
            self.qve.WriteOutHistograms()            
        
        elif self.mode == 'tds0-studies':
            outfilename=f'{self.outpath}/splitfiles/run{self.runNr}/tds0-studies.root'
            f=ROOT.TFile.Open(outfilename, 'recreate')
            for h in self.hists:
                hist=self.hists[h]
                f.WriteObject(hist, hist.GetName(), 'kOverwrite')
            print(f'Hists written to {outfilename}')
            f.Close()

    def LoadSystemObservables(self):
        systemobservablesfilename=f'{self.afswork}SystemObservables/run{self.TWCorrectionRun}/systemparameters.json'
        with open(systemobservablesfilename, 'r') as f:
            self.systemobservables=json.load(f)

    def yresidual3(self, detID):
        
        # self.MuFilter.GetPosition(detID,self.A,self.B)

        doca=self.muAna.Getyresidual(detID)
        s,p,b=self.muAna.parseDetID(detID)
        
        if s==3: return True
        key=10*s+p
        if key==12: return True # No hist for veto 3 yet

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

    def TimingDiscriminantCut(self):
        if self.timealignment=='old':
            if self.td<0: return False 
            else: return True
        else: 
            if not 'timingdiscriminant' in self.cutdists:
                print(f'No timing discriminant histogram')
                return 
            hist=self.cutdists['timingdiscriminant']
            mean, stddev=hist.GetMean(), hist.GetStdDev()
            if self.td < mean-2*stddev or self.td > mean+2*stddev: return False
            else: return True

    def dycut(self, hist, nsig=1):    
        dymin=hist.GetMean()-nsig*hist.GetStdDev()
        dymax=hist.GetMean()+nsig*hist.GetStdDev()
        return dymin, dymax

    def slopecut(self, slopecut=0.1):
        slopeX, slopeY = self.mom.x()/self.mom.z(), self.mom.y()/self.mom.z()
        if abs(slopeX)>slopecut or abs(slopeY)>slopecut: return 0
        else: return 1

    def tds0cut(self, nsig=1):
        if self.referencesystem==1: return 1
        
        if self.reft < 18 or self.reft > 21: return 0 
        else: return 1
