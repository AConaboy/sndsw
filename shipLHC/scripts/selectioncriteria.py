#!/usr/bin/env python

import ROOT
from pathlib import Path
from random import randint

A, B, locA, locB = ROOT.TVector3(), ROOT.TVector3(), ROOT.TVector3(), ROOT.TVector3()

class MuonSelectionCriteria(object):
    def Init(self, options, tw):
       
        self.tw=tw
        self.muAna=tw.muAna
        self.state=tw.state
        self.runNr=tw.r
        
        self.cuts={'OneHitPerSystem':options.OneHitPerSystem, 'SlopesCut':options.SlopesCut, 'nSiPMsCut':options.nSiPMsCut}

        if options.datalocation=='physics':
            self.outpath=f'{options.afswork}-physics2022/'
        elif options.datalocation=='commissioning':
            self.outpath=f'{options.afswork}-commissioning/'
        elif options.datalocation=='H8':
            self.outpath=f'{options.afswork}-H8/'

        self.EventNumber=-1

        self.systemAndPlanes = {1:2,2:5,3:4}
        self.systemAndBars = {1:7,2:10,3:60}
        self.nchs={1:224, 2:800}
        self.subsystemdict={1:'Veto', 2:'US', 3:'DS'}
        self.zPos=self.M.zPos

        self.freq=160.316E6
        self.TDC2ns=1E9/self.freq

        self.hists=self.M.h

        self.histtypes=['slopes', 'dy', 'nSiPMs']
        
        self.tracksinfo={}

        self.cutdists=self.muAna.GetCutDistributions(self.runNr, ('dy', 'timingdiscriminant'), task='SelectionCriteria')
        
        slopetitle = f'Track slopes;slope x [rad];slope y [rad]'
        self.hists[f'slopes'] = ROOT.TH2F(f'slopes',slopetitle, 300, -1.5, 1.5, 300, -1.5, 1.5)        
        
        slopetitle = f'Number of fired DS horizontal bars;Fired DS horizontal bars;Counts'
        self.hists[f'firedDSHbars'] = ROOT.TH1I(f'firedDSHbars',slopetitle,10, 0, 9)


        systemhistcriteria={None:None, '1DSb':'1 DS bar per plane', '1USb':'1 US bar per plane'}
        for s in (1, 2):
            for criterion in systemhistcriteria.keys():
                if not criterion:
                    title=self.subsystemdict[s]+' channel hit rate;Channel;Counts'
                    self.hists[f'channelhitrate-{self.subsystemdict[s]}']=ROOT.TH1F(f'channelhitrate-{self.subsystemdict[s]}', title, self.nchs[s]+1, 0, self.nchs[s])
                    
                    title='#splitline{'+self.subsystemdict[s]+' channel hit rate}{with timing discriminant cut};Channel;Counts'
                    self.hists[f'channelhitrate-{self.subsystemdict[s]}-tdcut']=ROOT.TH1F(f'channelhitrate-{self.subsystemdict[s]}-tdcut', title, self.nchs[s]+1, 0, self.nchs[s])

                    title='#splitline{'+self.subsystemdict[s]+' channel hit rate}{with dy cut};Channel;Counts'
                    self.hists[f'channelhitrate-{self.subsystemdict[s]}-dycut']=ROOT.TH1F(f'channelhitrate-{self.subsystemdict[s]}-dycut', title, self.nchs[s]+1, 0, self.nchs[s])

                    if s==2:
                        title='#chi_{#nu}^{2} of DS track fitting;#chi_{#nu}^{2} [dimensionless];Counts'
                        self.hists['trackfitting']=ROOT.TH1F('trackfitting', title, 50,0,10)

                        title='Prob. of higher #chi_{#nu}^{2} of DS track fitting;#chi_{#nu}^{2} [dimensionless];Counts'
                        self.hists['Probtrackfitting']=ROOT.TH1F('Probtrackfitting', title, 100,0,1)                        

                        timingdistitle="DS3H average - US1 average;t^{DS}_{0} - #bar{t_{US1}} [ns];Counts"
                        self.hists['timingdiscriminant']=ROOT.TH1F('timingdiscriminant',timingdistitle, 500,-25,50)        
                        
                        us1title="US1 average;#bar{t_{US1}} [ns];Counts"
                        self.hists['averageUS1time']=ROOT.TH1F('averageUS1time',us1title, 500,-25,50)              
                        
                        DS3Htitle="DS3H average time;t^{DS}_{0} [ns];Counts"
                        self.hists['averageDS3Htime']=ROOT.TH1F('averageDS3Htime',DS3Htitle, 500,-25,50)                    
                else: 
                    
                    title='#splitline{'+self.subsystemdict[s]+' channel hit rate}{with '+systemhistcriteria[criterion]+'};Channel;Counts'
                    self.hists[f'channelhitrate-{self.subsystemdict[s]}-{criterion}']=ROOT.TH1F(f'channelhitrate-{self.subsystemdict[s]}-{criterion}', title, self.nchs[s]+1, 0, self.nchs[s])
                
                    if s==2:   

                        title='#splitline{Prob. of higher #chi_{#nu}^{2} of DS track fitting}{with '+systemhistcriteria[criterion]+'};#chi_{#nu}^{2} [dimensionless];Counts'
                        self.hists[f'Probtrackfitting-{criterion}']=ROOT.TH1F(f'Probtrackfitting-{criterion}', title, 100,0,1)   

                        title='#splitline{#chi_{#nu}^{2} of DS track fitting}{with '+systemhistcriteria[criterion]+'};#chi_{#nu}^{2} [dimensionless];Counts'
                        self.hists[f'trackfitting-{criterion}']=ROOT.TH1F(f'trackfitting-{criterion}', title, 50,0,10)

                        timingdistitle="#splitline{DS3H average - US1 average}{with "+systemhistcriteria[criterion]+"};t^{DS}_{0} - #bar{t_{US1}} [ns];Counts"
                        self.hists[f'timingdiscriminant-{criterion}']=ROOT.TH1F(f'timingdiscriminant-{criterion}',timingdistitle, 500,-25,50)                    
                        
                        us1title="#splitline{US1 average}{with "+systemhistcriteria[criterion]+"};#bar{t_{US1}} [ns];Counts"
                        self.hists[f'averageUS1time-{criterion}']=ROOT.TH1F(f'averageUS1time-{criterion}',us1title, 500,-25,50)
                        
                        DS3Htitle="#splitline{DS3H average time}{with "+systemhistcriteria[criterion]+"};t^{DS}_{0} [ns];Counts"
                        self.hists[f'averageDS3Htime-{criterion}']=ROOT.TH1F(f'averageDS3Htime-{criterion}',DS3Htitle, 500,-25,50)
                        
        timingdisttitle="#splitline{DS3H average - US1 average}{using median for US1 left/right time};t^{DS}_{0} - #bar{t_{US1}} [ns];Counts"
        self.hists[f'timingdiscriminant-{criterion}']=ROOT.TH1F(f'timingdiscriminant-medianUS1left/right',timingdisttitle, 500,-25,50)                    

        NPlanestitle='Fraction of expected planes firing with only 1 scintillator;N_{fired}/N_{expected} [dimensionless];Counts'
        self.hists['NPlanes']=ROOT.TH1F('NPlanes', NPlanestitle, 10, 0., 1.)

        NRecoTrackstitle='Number of tracks reconstructed in Scifi and DS;N tracks [dimensionless];Counts'
        self.hists['NRecoTracks']=ROOT.TH1F('NRecoTracks', NRecoTrackstitle, 5, 0, 5)        

        for criterion in ('OneHitPerSystem', 'Atleast2planeswith1scintillator', None):
            chi2xypredtitle='#splitline{Fitted DS track #chi^{2}_{#nu}# correlated with track xy position}{evaluated at DS '+str(options.chi2xpred_zpos+1)+'};x [cm];y [cm];DS track #chi^{2}_{#nu}#, [dimensionless]'
            if criterion: name=f'chi2xypred-{criterion}'
            else: name='chi2xypred'
            self.hists[name]=ROOT.TH3F(name, chi2xypredtitle, 100, -90, 10, 80, 0, 80, 400, 0, 400)

        self.hists['TDS0']=ROOT.TH1F('TDS0','Average of fired DS horizontal SiPMs;DS horizontal average [ns];Counts', 200, 0, 50)

        for i in (3,2):
            name=f'delta{i}{i-1}'
            title='Time difference between horizontal SiPMs in DS horizontal plane '+str(i)+' and '+str(i-1)+';DS'+str(i)+' - DS'+str(i-1)+' [ns]; Counts'
            self.hists[name]=ROOT.TH1F(name, title, 80, -10, 10)

        for subsystem in (1,2):
            for plane in range(self.systemAndPlanes[subsystem]):
                
                yrestitle=f'y-residual_{10*subsystem+plane};track y - bar y-midpoint [cm];track y-residual;Counts / 1 cm'
                self.hists[f'dy_{10*subsystem+plane}']=ROOT.TH1F(f'dy_{10*subsystem+plane}', yrestitle, 60, -30, 30)

                nSiPMstitle=f'nSiPMs_{10*subsystem+plane};nSiPMs left;nSiPMs right'
                self.hists[f'nSiPMs_{10*subsystem+plane}']=ROOT.TH2F(f'nSiPMs_{10*subsystem+plane}', nSiPMstitle,9,0,9,9,0,9)
        
        if not self.cuts['OneHitPerSystem']: # Make hist for correlating plane multiplicity and bar number
            for subsystem in (1,2):
                    title=f'Correlation between {self.subsystemdict[subsystem]} plane hit multiplicity and bar number;Plane multiplicity [dimensionless];Bar number [dimensionless]'
                    self.hists[f'{subsystem}_multiplicity']=ROOT.TH2I(f'{subsystem}_multiplicity', title, 10,0,10,10,0,10)           
        if not self.cuts['OneHitPerSystem']:
            for subsystem in (1,2):
                for plane in range(self.systemAndPlanes[subsystem]):
                    title=f'{self.subsystemdict[subsystem]} plane {plane} correlation between fired bars;Bar number, bar closest to track; Additional bar number'
                    self.hists[f'{10*subsystem+plane}_barcorrelation']=ROOT.TH2I(f'{10*subsystem+plane}_barcorrelation', title, 10,0,9,10,0,9)
                    self.hists[f'{10*subsystem+plane}_barcorrelation'].GetXaxis().SetTitleOffset(0.8)
                    self.hists[f'{10*subsystem+plane}_barcorrelation'].GetXaxis().SetTitleSize(0.05)
                    self.hists[f'{10*subsystem+plane}_barcorrelation'].GetYaxis().SetTitleSize(0.05)

        for plane in range(7):
            if plane in (0,2,4): pull, coord, station='Vertical', 'y', int(plane/2)
            elif plane in (1,3,5): pull, coord, station='Horizontal', 'x', int((plane-1)/2)
            else: pull, coord, station='Horizontal', 'x', 3
            title=f'{pull} pull for DS plane {station};Bar position in {coord} [cm];Expected position in {coord} [cm]'
            if coord=='y': bins=(80, 0, 80, 80, 0, 80)
            else: bins=(70, -60, 10, 70, -60, 10)
            self.hists[f'pull_{3*10+plane}']=ROOT.TH2F(f'pull_{plane}', title, *bins)

        for plane in range(self.systemAndPlanes[2]):
            for bar in range(self.systemAndBars[2]):
                title='DS track x_{predicted} against Scifi track x_{predicted} for tracks extrapolated to US plane '+str(plane+1)+', bar '+str(bar+1)+';DS track x_{predicted} - Scifi track x_{predicted}[cm]'
                # self.hists[f'DSxvScifix_plane{plane}_bar{bar}']=ROOT.TH2F(f'DSxvScifix_plane{plane}_bar{bar}', title,  100, -90, 10, 100, -90, 10)
                self.hists[f'DSxvScifix_plane{plane}_bar{bar}']=ROOT.TH1F(f'DSxvScifix_plane{plane}_bar{bar}', title,  40, -20, 20)
        
        for subsystem in (1,2):
            for plane in range(self.systemAndPlanes[subsystem]):
                if subsystem==3 and plane==3: continue # No 4th horizontal plane
                for bar in range(self.systemAndBars[subsystem]):
                    title=self.subsystemdict[subsystem]+' plane '+str(plane+1)+' bar '+str(bar+1)+' average of median '+self.state+' time from each side as a function of x_{predicted};x_{predicted} [cm];'+'#frac{1}{2}#times(t_{left}+t_{right}) [ns];Counts'
                    self.hists[f'{10*subsystem+plane}_bar{bar}_averagetime_{self.state}']=ROOT.TH2F(f'{10*subsystem+plane}_bar{bar}_averagetime_{self.state}', title, 100,0,100,250,0,25)
                    self.hists[f'{10*subsystem+plane}_bar{bar}_averagetime_{self.state}'].GetXaxis().SetTitleOffset(0.8)
                    self.hists[f'{10*subsystem+plane}_bar{bar}_averagetime_{self.state}'].GetXaxis().SetTitleSize(0.05)
                    self.hists[f'{10*subsystem+plane}_bar{bar}_averagetime_{self.state}'].GetYaxis().SetTitleSize(0.05)

    def ExecuteEvent(self, event):
        
        hists=self.hists
        hits=event.Digi_MuFilterHits

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
        self.fstate=dstrack.getFittedState()
        
        self.pos=self.fstate.getPos()
        self.mom=self.fstate.getMom()
        self.trackchi2NDF=fitStatus.getChi2()/fitStatus.getNdf()+1E-10

        self.OneHitPerUS=self.muAna.OneHitPerSystem(hits, (2,))
        self.OneHitPerDS=self.muAna.OneHitPerSystem(hits, (3,))        

        if 'slopes' in self.histtypes: 
            self.slopes=self.mom.x()/self.mom.z(), self.mom.y()/self.mom.z()
            hists[f'slopes'].Fill(*self.slopes)

        # If fitted DS track passes slope criteria. Then fill correlate red-chi2 with xy-position
        if self.cuts['SlopesCut'] and not self.slopecut(): return
        
        TDS0, firedDSHbars=self.muAna.GetDSHaverage(self.MuFilter, hits) # Now returns in ns
        if TDS0==-999: print(f'Event {self.M.EventNumber} has a fitted DS track but no DS horizontal bars with both SiPMs firing.')
        self.hists['TDS0'].Fill(TDS0)
        self.hists['firedDSHbars'].Fill(firedDSHbars)

        # Use the correct subsystems for 1 bar / plane cut
        if inVeto: tmp=(1,2,3)
        else: tmp=(2,3)
        Nfired=self.muAna.OneHitPerSystem(hits, tmp, Nfired=True)
        hists['NPlanes'].Fill(Nfired)        

        self.redchi2hists(fitStatus)        

        if self.options.debug:
            self.FillChannelRateHists()
    
            self.Fillchi2xy(tmp)

        # if self.cuts['OneHitPerSystem']:
        self.FillTimingDiscriminantHists(hits, tmp)
        
        # Both of these methods loop through the hits
        self.FillPlaneMultiplicity()
        self.FillFiredBarCorrelation()

        # Loops through hits
        self.FillPullHists()
        
        # Loops through hits
        self.FillnSiPMs()
        
        # Loops through hits
        self.Fillyresidual()
        
        self.FilldeltaDSH()
        
        # Loops through hits, requires at least 6 SiPMs per side in Veto and 4 large SiPMs per side in US
        self.FillAverageBarTime()

        if inVeto and inDS:
            self.FillScifiDSresidual()
             
        if randint(0,1) and self.options.WriteOutTrackInfo:
            self.WriteOutTrack()            

    # def Fillyresidual(self,detID):
    def Fillyresidual(self):
    
        hits=self.M.eventTree.Digi_MuFilterHits
        for hit in hits:
            detID=hit.GetDetectorID()
            s,p,b=self.muAna.parseDetID(detID)
            if s==3: continue
            
            # zEx=self.zPos['MuFilter'][s*10+p]
            # lam=(zEx-self.pos.z())/self.mom.z()
            # self.Ex=ROOT.TVector3(self.pos.x()+lam*self.mom.x(), self.pos.y()+lam*self.mom.y(), self.pos.z()+lam*self.mom.z())
            self.MuFilter.GetPosition(detID, A,B)
            
            # MuonSelectionCriteria.Getyresidual(detID) requires that A has been filled already with the coordinates of detID! 
            pq = A-self.pos
            uCrossv= (B-A).Cross(self.mom)
            doca = pq.Dot(uCrossv)/uCrossv.Mag()
            self.hists[f'dy_{10*s+p}'].Fill(doca)
    
    def Getyresidual(self, detID):
        pq = A-self.pos
        uCrossv= (B-A).Cross(self.mom)
        doca = pq.Dot(uCrossv)/uCrossv.Mag()
        return doca   

    def yresidual3(self, detID):
        self.MuFilter.GetPosition(detID, A, B)
        doca=self.Getyresidual(detID)
        s,p,b=self.muAna.parseDetID(detID)
        if s==3: return True
        key=10*s+p
        
        dy_min, dy_max = MuonSelectionCriteria.dycut(self.cutdists[f'dy_{key}'])
        if doca>dy_max or doca<dy_min: return False
        else: return True

    @staticmethod
    def dycut(hist, nsig=1):    
        dymin=hist.GetMean()-nsig*hist.GetStdDev()
        dymax=hist.GetMean()+nsig*hist.GetStdDev()
        return dymin, dymax              

    def Fillchi2xy(self, tmp, key=30):
        z=self.zPos['MuFilter'][key]
        lam=(z-self.pos.z())/self.mom.z()
        Ex=ROOT.TVector3(self.pos.x()+lam*self.mom.x(), self.pos.y()+lam*self.mom.y(), self.pos.z()+lam*self.mom.z())
        x, y = Ex.x(), Ex.y()
        
    def slopecut(self, slopecut=0.1):
        slopeX, slopeY = self.slopes
        if abs(slopeX)>slopecut or abs(slopeY)>slopecut: return 0
        else: return 1

    def FillnSiPMs(self):
        hits=self.M.eventTree.Digi_MuFilterHits
        for hit in hits:
            detID=hit.GetDetectorID()
            nLeft, nRight=self.muAna.GetnFiredSiPMs(hit)
            s,p,b=self.muAna.parseDetID(detID)
            if s==3: continue
            key=10*s+p
            self.hists[f'nSiPMs_{key}'].Fill(nLeft, nRight)

    def nSiPMscut(self, hit, nLeft, nRight):
        s,p,b=self.muAna.parseDetID(hit.GetDetectorID())
        if s==1: 
            if nLeft<6 or nRight<6: return False
        elif s==2:
            if nLeft<4 or nRight<4: return False
        elif s==3: 
            pass 
        return True

    def FillPullHists(self):
        hits=self.M.eventTree.Digi_MuFilterHits
        for hit in hits:
            vertical=False
            detID=hit.GetDetectorID()
            s,p,b=self.muAna.parseDetID(detID)
            if s!=3: continue
            # parseDetID returns station number, zPos distinguishes between the planes in DS stations. 
            # Analysis.GetDSPlaneNumber(detID) returns the plane number to use in the zPos dictionary
            plane=self.muAna.GetDSPlaneNumber(detID)
            key=10*s+plane
            z=self.zPos['MuFilter'][key]
            self.MuFilter.GetPosition(detID, A, B)
            if self.muAna.DSVcheck(detID): vertical=True 
            
            if vertical: 
                barposition=0.5*(A.x()+B.x()) 
                lam=(z-self.pos.z())/self.mom.z()
                tmp=ROOT.TVector3(self.pos.x()+lam*self.mom.x(), self.pos.y()+lam*self.mom.y(), self.pos.z()+lam*self.mom.z())             
                histname=f'pull_{key}'
                self.hists[histname].Fill(barposition, tmp.x())  
            
            else: 
                barposition=0.5*(A.y()+B.y()) 
                lam=(z-self.pos.z())/self.mom.z()
                tmp=ROOT.TVector3(self.pos.x()+lam*self.mom.x(), self.pos.y()+lam*self.mom.y(), self.pos.z()+lam*self.mom.z())             
                histname=f'pull_{key}'
                self.hists[histname].Fill(barposition, tmp.y())

    def redchi2hists(self, fitStatus):

        self.hists[f'trackfitting'].Fill(self.trackchi2NDF)
        chi2, ndf=fitStatus.getChi2(), fitStatus.getNdf()
        # print(f'chi2: {chi2}, prob: {ndf}')
        if ndf!=0: 
            prob=ROOT.TMath.Prob(chi2, int(ndf))
            self.hists[f'Probtrackfitting'].Fill(prob)       
        if self.OneHitPerDS: 
            self.hists[f'trackfitting-1DSb'].Fill(self.trackchi2NDF)       
            if ndf!=0:self.hists[f'Probtrackfitting-1DSb'].Fill(prob)       
        if self.OneHitPerUS: 
            self.hists[f'trackfitting-1USb'].Fill(self.trackchi2NDF)
            if ndf!=0:self.hists[f'Probtrackfitting-1USb'].Fill(prob)
         
    def FillAverageBarTime(self, state=None):
        if not state: state=self.state
        
        hits=self.M.eventTree.Digi_MuFilterHits
        for hit in hits:
            detID=hit.GetDetectorID()
            s,p,b=self.muAna.parseDetID(detID)
            if s==3: continue
            medians=self.muAna.GetMedianTime(hit, mode=state)
            if not medians: continue
            if len(medians)!=2: continue
            
            zEx=self.zPos['MuFilter'][s*10+p]
            lam=(zEx-self.pos.z())/self.mom.z()
            self.Ex=ROOT.TVector3(self.pos.x()+lam*self.mom.x(), self.pos.y()+lam*self.mom.y(), self.pos.z()+lam*self.mom.z())
            self.MuFilter.GetPosition(detID, A, B)
            xpred=ROOT.TMath.Sqrt((A.x()-self.Ex.x())**2+(A.y()-self.Ex.y())**2+(A.z()-self.Ex.z())**2)
            
            averagetime=sum(medians.values())/2 
            hist=self.hists[f'{s*10+p}_bar{b}_averagetime_{self.state}']
            hist.Fill(xpred, averagetime)
            
    def FillPlaneMultiplicity(self):
        hits=self.M.eventTree.Digi_MuFilterHits
        multiplicity = {s:{i:[] for i in range(self.systemAndPlanes[s])} for s in (1,2)}

        for hit in hits:
            detID=hit.GetDetectorID()
            s,p,b = self.muAna.parseDetID(detID)
            if s==3: continue
            multiplicity[s][p].append(b)
        for s in multiplicity:
            for p in multiplicity[s]:
                bs=multiplicity[s][p]
                for b in bs: self.hists[f'{s}_multiplicity'].Fill(len(bs), b)
                
                
    def FillTimingDiscriminantHists(self, hits, tmp):
        
        timingdiscriminant=self.GetTimingDiscriminant()
        if timingdiscriminant==-420:
            return
        
        self.hists['timingdiscriminant'].Fill(timingdiscriminant)
        
        timingdistfrommedian=self.GetTimingDiscriminant(mode='median')

        ds3haverage=self.muAna.GetDSHaverage(self.MuFilter, hits, mode='timingdiscriminant')
        if ds3haverage>-420: 
            us1averagetime = -1*(timingdiscriminant - ds3haverage)
            self.hists['averageDS3Htime'].Fill(ds3haverage)
            self.hists['averageUS1time'].Fill(us1averagetime)
            if self.OneHitPerUS: 
                self.hists['timingdiscriminant-1USb'].Fill(timingdiscriminant)
                self.hists['averageDS3Htime-1USb'].Fill(ds3haverage)
                self.hists['averageUS1time-1USb'].Fill(us1averagetime)
            if self.OneHitPerDS:
                self.hists['timingdiscriminant-1DSb'].Fill(timingdiscriminant)
                self.hists['averageDS3Htime-1DSb'].Fill(ds3haverage)
                self.hists['averageUS1time-1DSb'].Fill(us1averagetime)
                    
                
    def GetTimingDiscriminant(self, mode='mean'):
        
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

        if mode=='mean':averageUS1time=self.muAna.GetAverageTime(us1hit)
        if mode=='median':
            medians=self.muAna.GetMedianTime(us1hit)
            if not medians: return -999
            averageUS1time=sum(medians.values())/len(medians)
        if not averageUS1time: return -998

        DS3Haverage=self.muAna.GetDSHaverage(self.MuFilter, hits, mode='timingdiscriminant')

        return DS3Haverage-averageUS1time                    
                
    def FillFiredBarCorrelation(self):
        hits=self.M.eventTree.Digi_MuFilterHits
        multiplicity = {s:{i:[] for i in range(self.systemAndPlanes[s])} for s in (1,2)}
        for hit in hits:
            subsystem, plane, bar = self.muAna.parseDetID(hit.GetDetectorID())
            if subsystem==3: continue
            multiplicity[subsystem][plane].append(bar)
        
        for subsystem in multiplicity:
            for plane in multiplicity[subsystem]:
                if len(multiplicity[subsystem][plane])<2: continue
                bars=multiplicity[subsystem][plane]
                docas={}
                for bar in bars:
                    detID=self.muAna.MakeDetID((subsystem, plane, bar))
                    self.MuFilter.GetPosition(detID, A, B)
                    docas[self.Getyresidual(detID)]=bar

                closest_bar=docas.pop(min(docas))
                for x in docas: self.hists[f'{10*subsystem+plane}_barcorrelation'].Fill(closest_bar, docas[x])

    def FilldeltaDSH(self):
        res=self.muAna.GetDSHaverage(self.MuFilter, self.M.eventTree.Digi_MuFilterHits, mode='deltastations')
        if not res: return
        [self.hists[i].Fill(res[i]) for i in res]

    def FillScifiDSresidual(self):
        scifi_track=self.M.Reco_MuonTracks[1]

        scifi_fstate=scifi_track.getFittedState()
        scifi_fitStatus=scifi_track.getFitStatus()
        if not scifi_fitStatus.isFitConverged(): return
        scifi_pos=scifi_fstate.getPos()
        scifi_mom=self.fstate.getMom()
        
        hits=self.M.eventTree.Digi_MuFilterHits
        for hit in hits:
            detID=hit.GetDetectorID()
            s,p,b=self.muAna.parseDetID(detID)
            if s!=2: continue
            
            zEx=self.zPos['MuFilter'][s*10+p]
            ds_lam=(zEx-self.pos.z())/self.mom.z()
            self.Ex=ROOT.TVector3(self.pos.x()+ds_lam*self.mom.x(), self.pos.y()+ds_lam*self.mom.y(), self.pos.z()+ds_lam*self.mom.z())
            scifi_lam=(zEx-scifi_pos.z())/scifi_mom.z()
            scifi_Ex=ROOT.TVector3(scifi_pos.x()+scifi_lam*scifi_mom.x(), scifi_pos.y()+scifi_lam*scifi_mom.y(), scifi_pos.z()+scifi_lam*scifi_mom.z())
            # self.MuFilter.GetPosition(detID, A, B)

            self.hists[f'DSxvScifix_plane{p}_bar{b}'].Fill(self.Ex.x() - scifi_Ex.x())
            
    def FillChannelRateHists(self):
        hits=self.M.eventTree.Digi_MuFilterHits
        for hit in hits:
            detID=hit.GetDetectorID()
            s,p,b=self.muAna.parseDetID(detID)

            if s==3: continue
            for x in hit.GetAllTimes():
                SiPM, clock = x
                channel=self.muAna.GetSiPMNumberInSystem_PCBbyPCB(detID, SiPM)
                if self.OneHitPerUS: self.hists[f'channelhitrate-{self.subsystemdict[s]}-1USb'].Fill(channel)
                if self.OneHitPerDS: self.hists[f'channelhitrate-{self.subsystemdict[s]}-1DSb'].Fill(channel)
                self.hists[f'channelhitrate-{self.subsystemdict[s]}'].Fill(channel)
                td=self.GetTimingDiscriminant()
                if not self.TimingDiscriminantCut(td): self.hists[f'channelhitrate-{self.subsystemdict[s]}-tdcut'].Fill(channel)
                if not self.yresidual3(detID): self.hists[f'channelhitrate-{self.subsystemdict[s]}-dycut'].Fill(channel)

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

    def WriteOutHistograms(self):

        outpath=f'{self.outpath}/splitfiles/run{self.runNr}/SelectionCriteria/'
        path_obj=Path(outpath)
        path_obj.mkdir(parents=True, exist_ok=True)        
        outfile=f'SelectionCriteria_{self.options.nStart}.root'
        f=ROOT.TFile.Open(outpath+outfile, 'update')
        
        additionalkeys=['averagetime', 'DSxvScifix']
        for histname in self.M.h:
            if not any( [histname.find(additionalkey)>-1 for additionalkey in additionalkeys] ):
                
                hist=self.M.h[histname]
                f.WriteObject(hist, hist.GetName(), 'kOverwrite')
            else:
                for additionalkey in additionalkeys:
                    if histname.find(additionalkey)>-1:
                        if not hasattr(f, additionalkey): residualsfolder=f.mkdir(additionalkey)
                        else: residualsfolder=f.Get(additionalkey)
                        residualsfolder.cd()
                        hist=self.M.h[histname]
                        hist.Write(histname, 2) # The 2 means it will overwrite a hist of the same name
        f.Close()
        print(f'{len(self.M.h)} histograms saved to {outpath}{outfile}')
        
    def WriteOutToTestFile(self):       
        outfile=f'testing_{self.options.runNumber}_{self.options.nStart}_{self.options.nEvents}.root'
        f=ROOT.TFile.Open(outfile, 'update')
        
        # mainhistkeys=['dy', 'slopes', 'nSiPMs', 'timingdiscriminant', 'trackfitting', 'TDS0', ]
        additionalhistkeys=['averagetime', 'DSxvScifix']
        for histname in self.M.h:
            if not any( [histname.find(additionalhistkey)>-1 for additionalhistkey in additionalhistkeys] ):
                hist=self.M.h[histname]
                f.WriteObject(hist, hist.GetName(), 'kOverwrite')
            else:
                for additionalkey in additionalhistkeys:
                    if histname.find(additionalkey)>-1:
                        if not hasattr(f, additionalkey): residualsfolder=f.mkdir(additionalkey)
                        else: residualsfolder=f.Get(additionalkey)
                        residualsfolder.cd()
                        hist=self.M.h[histname]
                        hist.Write(histname, 2)

        f.Close()
        print(f'{len(self.M.h)} histograms saved to {outfile}')   
    
    def WriteOutTrack(self):
        
        tracks=self.M.Reco_MuonTracks 
        if len(tracks)==1:  DStrack=self.M.Reco_MuonTracks[0]
        else:
            for idx, tr in enumerate(self.M.Reco_MuonTracks):
                if tr.GetUniqueID()==3:
                    DStrack=self.M.Reco_MuonTracks[idx]
        
        self.fstate=DStrack.getFittedState()
        mom, pos=self.fstate.getMom(), self.fstate.getPos()
        detIDs=[hit.GetDetectorID() for hit in self.M.eventTree.Digi_MuFilterHits]
        
        trackinfo={'mom':[mom.x(), mom.y(), mom.z()], 'pos':[pos.x(),pos.y(),pos.z(),], 'detIDs':detIDs}
        self.tracksinfo[self.M.EventNumber]=trackinfo
        
    def SaveTrackInfos(self):
        import json 
        with open(f'{self.outpath}Results/run{self.runNr}/trackinfos.json', 'w') as outfile:
            json.dump(self.tracksinfo, outfile)