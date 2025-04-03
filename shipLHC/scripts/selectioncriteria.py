#!/usr/bin/env python

import ROOT
from pathlib import Path
from random import randint

A, B, locA, locB = ROOT.TVector3(), ROOT.TVector3(), ROOT.TVector3(), ROOT.TVector3()

class MuonSelectionCriteria(object):
    def __init__(self, options, tw):
       
        self.tw=tw
        self.options=options
        self.muAna=tw.muAna
        self.state=tw.state
        self.runNr=tw.runNr
        
        self.cuts={'OneHitPerSystem':options.OneHitPerSystem, 'SlopesCut':options.SlopesCut}

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
        self.zPos=tw.zPos

        self.freq=160.316E6
        self.TDC2ns=1E9/self.freq

        self.MuFilter=tw.MuFilter
        self.Scifi=tw.Scifi

        self.hists=tw.hists

        self.histtypes=['slopes', 'dy']
        
        self.tracksinfo={}

        # self.cutdists=self.muAna.GetCutDistributions(self.runNr, ('dy', 'timingdiscriminant'), task='SelectionCriteria')
        
        slopetitle = f'Track slopes;slope x [rad];slope y [rad]'
        self.hists[f'slopes'] = ROOT.TH2F(f'slopes',slopetitle, 300, -1.5, 1.5, 300, -1.5, 1.5)        

        slopetitle = f'Number of fired DS horizontal bars;Fired DS horizontal bars;Counts'
        self.hists[f'firedDSHbars'] = ROOT.TH1I(f'firedDSHbars',slopetitle,10, 0, 9)

        for p in ['x', 'y']:
            spatialchi2pNDF = '#chi^{2}_{#nu} of fitted tracks\ncorrelated with '+p+' position at the first muon system plane;#chi^{2}_{#nu} [dimensionless];'+p+' [cm]'
            spatial_binning = (110, -100, 10) if p=='x' else (80,0,80)
            self.hists[f'chi2pNDFv{p}'] = ROOT.TH2F(f'chi2pNDFv{p}', spatialchi2pNDF, 20, 0, 10, *spatial_binning)

        cutflowtitle = f'Cut flow;Cut name;Events that pass cut'
        
        # cuts = ['has track', '#chi^{2}_{#nu}', 't^{ds}_{0}', 'slope']
        cuts = ['has track', '#chi^{2}_{#nu}', 'slope']
        
        self.cuteffect = {cut: 0 for cut in cuts}
        self.tds0vcuts = {cut: ROOT.TH1D(f'{cut}-tds0', 't_{0}^{DS} v '+cut+' cut', 100, 0, 25) for cut in cuts}
        self.deltads32vcuts = {cut: ROOT.TH1D(f'{cut} #Delta ds32', '#Delta (t_{3}^{DS}, t_{2}^{DS}) v '+cut+' cut', 80, -10, 10) for cut in cuts}
        self.deltads21vcuts = {cut: ROOT.TH1D(f'{cut} #Delta ds21', '#Delta (t_{2}^{DS}, t_{1}^{DS}) v '+cut+' cut', 80, -10, 10) for cut in cuts}

        self.hists[f'cutflow'] = ROOT.TH1I(f'cutflow',cutflowtitle,len(self.cuteffect), 0, len(self.cuteffect))

        s=2
        systemhistcriteria={None:None, 'slopecut':'slope cut'}
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

        for i in (3,2):
            name=f'delta{i}{i-1}'
            title='Time difference between horizontal SiPMs in DS horizontal plane '+str(i)+' and '+str(i-1)+';t(DS'+str(i)+') - t(DS'+str(i-1)+') [ns]; Counts'
            self.hists[name]=ROOT.TH1F(name, title, 80, -10, 10)

        for st in range(5):
            name=f'scifi-delta{st+1}{st}'
            title='Time difference between average time measured in Scifi plane '+str(st+1)+' and '+str(st)+';#Delta(station' +str(st+1)+','+str(st)+') [ns];Counts'
            self.hists[name]=ROOT.TH1F(name, title, 80, -10, 10)

        for plane in range(self.systemAndPlanes[s]):
            
            yrestitle=f'y-residual_{10*s+plane};track y - bar y-midpoint [cm];Counts / 1 cm'
            self.hists[f'dy_{10*s+plane}']=ROOT.TH1F(f'dy_{10*s+plane}', yrestitle, 60, -30, 30)

            nSiPMstitle=f'nSiPMs_{10*s+plane};nSiPMs left;nSiPMs right'
            self.hists[f'nSiPMs_{10*s+plane}']=ROOT.TH2F(f'nSiPMs_{10*s+plane}', nSiPMstitle,9,0,9,9,0,9)
        
        if not self.cuts['OneHitPerSystem']: # Make hist for correlating plane multiplicity and bar number
            for subsystem in (1,2):
                    title=f'Correlation between {self.subsystemdict[s]} plane hit multiplicity and bar number;Plane multiplicity [dimensionless];Bar number [dimensionless]'
                    self.hists[f'{s}_multiplicity']=ROOT.TH2I(f'{s}_multiplicity', title, 10,0,10,10,0,10)           
        if not self.cuts['OneHitPerSystem']:
            for plane in range(self.systemAndPlanes[s]):
                title=f'{self.subsystemdict[s]} plane {plane} correlation between fired bars;Bar number, bar closest to track; Additional bar number'
                self.hists[f'{10*s+plane}_barcorrelation']=ROOT.TH2I(f'{10*s+plane}_barcorrelation', title, 10,0,9,10,0,9)
                self.hists[f'{10*s+plane}_barcorrelation'].GetXaxis().SetTitleOffset(0.8)
                self.hists[f'{10*s+plane}_barcorrelation'].GetXaxis().SetTitleSize(0.05)
                self.hists[f'{10*s+plane}_barcorrelation'].GetYaxis().SetTitleSize(0.05)

        for plane in range(7):
            if plane in (0,2,4): pull, coord, station='Vertical', 'y', int(plane/2)
            elif plane in (1,3,5): pull, coord, station='Horizontal', 'x', int((plane-1)/2)
            else: pull, coord, station='Horizontal', 'x', 3
            title=f'{pull} pull for DS plane {station};Bar position - expected position in {coord} [cm];Number of fired bars in muon system cluster'
            # if coord=='y': bins=(80, 0, 80, 80, 0, 80)
            # else: bins=(70, -60, 10, 70, -60, 10)
            bins=(30, -15, 15), (6,0,6)
            self.hists[f'pull_{3*10+plane}']=ROOT.TH2F(f'pullvclustermultiplicity_{plane}', title, 30, -15, 15, 6, 0, 6)

            # Make hists of the times in each DS plane
            name=f'DSplane{plane}-time'
            title=f'Times measured in DS plane {plane+1};Time [ns];Counts'
            self.hists[name]=ROOT.TH1F(name, title, 100, 0, 25)

        for plane in range(self.systemAndPlanes[2]):
            for bar in range(self.systemAndBars[2]):
                title='DS track x_{predicted} against Scifi track x_{predicted} for tracks extrapolated to US plane '+str(plane+1)+', bar '+str(bar+1)+';DS track x_{predicted} - Scifi track x_{predicted}[cm]'
                # self.hists[f'DSxvScifix_plane{plane}_bar{bar}']=ROOT.TH2F(f'DSxvScifix_plane{plane}_bar{bar}', title,  100, -90, 10, 100, -90, 10)
                self.hists[f'DSxvScifix_plane{plane}_bar{bar}']=ROOT.TH1F(f'DSxvScifix_plane{plane}_bar{bar}', title,  40, -20, 20)
        
        # for plane in range(self.systemAndPlanes[s]):
        #     if subsystem==3 and plane==3: continue # No 4th horizontal plane
        #     for bar in range(self.systemAndBars[s]):
        #         title=self.subsystemdict[s]+' plane '+str(plane+1)+' bar '+str(bar+1)+' average of median '+self.state+' time from each side as a function of x_{predicted};x_{predicted} [cm];'+'#frac{1}{2}#times(t_{left}+t_{right}) [ns];Counts'
        #         self.hists[f'{10*s+plane}_bar{bar}_averagetime_{self.state}']=ROOT.TH2F(f'{10*subsystem+plane}_bar{bar}_averagetime_{self.state}', title, 100,0,100,250,0,25)
        #         self.hists[f'{10*s+plane}_bar{bar}_averagetime_{self.state}'].GetXaxis().SetTitleOffset(0.8)
        #         self.hists[f'{10*s+plane}_bar{bar}_averagetime_{self.state}'].GetXaxis().SetTitleSize(0.05)
        #         self.hists[f'{10*s+plane}_bar{bar}_averagetime_{self.state}'].GetYaxis().SetTitleSize(0.05)
    
        # self.reft, firedDSHbars=self.muAna.GetDSHaverage(self.MuFilter, hits) # Now returns in ns
        # if self.reft==-999: print(f'Event {self.M.EventNumber} has a fitted DS track but no DS horizontal bars with both SiPMs firing.')
        # self.hists['reft'].Fill(self.reft)
        # self.hists['firedDSHbars'].Fill(firedDSHbars)      

    def FillHists(self, hits, scifi_hits):
        
        # Loops through hits
        self.FillnSiPMs(hits)

        if self.options.debug and self.tw.hasTrack:
            self.FillChannelRateHists()
    
        # if self.cuts['OneHitPerSystem']:
        if self.tw.hasTrack and self.tw.referencesystem==3:
            self.FillTimingDiscriminantHists(hits)
            self.Fillyresidual(hits)
            self.FilldeltaDSH(hits)   
            self.FillSlopeHists(hits)
            self.FillPullHists(hits)
            self.Fillchi2xy()
            # self.FillScifiDSresidual()

        elif self.tw.hasTrack and self.tw.referencesystem==1:
            self.FilldeltaSF(scifi_hits)
        
        # Both of these methods loop through the hits
        # self.FillPlaneMultiplicity(hits)
        # self.FillFiredBarCorrelation(hits)
   
        # Loops through hits, requires at least 6 SiPMs per side in Veto and 4 large SiPMs per side in US
        # self.FillAverageBarTime(hits)
             
        # if randint(0,1) and self.options.WriteOutTrackInfo:
        #     self.WriteOutTrack()            

    def Fillyresidual(self, hits):
    
        for hit in hits:
            detID=hit.GetDetectorID()
            s,p,b=self.muAna.parseDetID(detID)
            if s!=2: continue
            
            self.tw.MuFilter.GetPosition(detID, A,B)
            
            # MuonSelectionCriteria.Getyresidual(detID) requires that A has been filled already with the coordinates of detID! 
            pq = A-self.tw.pos
            uCrossv= (B-A).Cross(self.tw.mom)
            doca = pq.Dot(uCrossv)/uCrossv.Mag()
            self.hists[f'dy_{10*s+p}'].Fill(doca)
    
    def Getyresidual(self, detID):
        pq = A-self.tw.pos
        uCrossv= (B-A).Cross(self.tw.mom)
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

    def FillSlopeHists(self, hits):
        self.slopes=self.tw.mom.x()/self.tw.mom.z(), self.tw.mom.y()/self.tw.mom.z()
        self.hists[f'slopes'].Fill(*self.slopes)

    @staticmethod
    def dycut(hist, nsig=1):    
        dymin=hist.GetMean()-nsig*hist.GetStdDev()
        dymax=hist.GetMean()+nsig*hist.GetStdDev()
        return dymin, dymax              

    def Fillchi2xy(self, key=30):
        z=self.zPos['MuFilter'][key]
        lam=(z-self.tw.pos.z())/self.tw.mom.z()
        Ex=ROOT.TVector3(self.tw.pos.x()+lam*self.tw.mom.x(), self.tw.pos.y()+lam*self.tw.mom.y(), self.tw.pos.z()+lam*self.tw.mom.z())
        x, y = Ex.x(), Ex.y()

        self.hists[f'chi2pNDFvx'].Fill(self.tw.trackchi2NDF, x)
        self.hists[f'chi2pNDFvy'].Fill(self.tw.trackchi2NDF, y)
        
    def slopecut(self, slopecut=0.1):
        slopeX, slopeY = self.slopes
        if abs(slopeX)>slopecut or abs(slopeY)>slopecut: return 0
        else: return 1

    def FillnSiPMs(self, hits):

        for hit in hits:
            detID=hit.GetDetectorID()
            nLeft, nRight=self.muAna.GetnFiredSiPMs(hit)
            s,p,b=self.muAna.parseDetID(detID)
            if s!=2: continue
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

    def FillPullHists(self, hits):
        ### Need to update this to only use hits used in the track fitting! Can copy code from the tds0 determination

        if not all([self.tw.passslopecut, self.tw.passredchi2cut]): return 

        stations = {i.GetDetectorID():i for i in hits if i.GetDetectorID()//10000==3}
        station_times = {station:[] for station in range(7)}
        # for hit in hits:
        track = self.tw.M.Reco_MuonTracks[0]
        nM = track.getNumPointsWithMeasurement()
        for n in range(nM):
            M=track.getPointWithMeasurement(n)
            W=M.getRawMeasurement()
            nbars = int(W.getRawHitCoords()[6])
            detID = W.getDetId()
            s,p,b=self.muAna.parseDetID(detID)
            if s!=3: continue

            # parseDetID returns station number, zPos distinguishes between the planes in DS stations. 
            # Analysis.GetDSPlaneNumber(detID) returns the plane number to use in the zPos dictionary

            plane=self.muAna.GetDSPlaneNumber(detID)
            key=10*s+plane
            z=self.zPos['MuFilter'][key]
            self.tw.MuFilter.GetPosition(detID, A, B)

            # Stores times for planes
            hit = stations[detID]
            # station_times[plane].append([i*self.TDC2ns for i in hit.GetAllTimes()])
            for x in hit.GetAllTimes():
                SiPM, cc = x
                dscorrectedtime=self.tw.MuFilter.GetCorrectedTime(detID, SiPM, cc*self.TDC2ns, 0)
                self.hists[f'DSplane{plane}-time'].Fill(dscorrectedtime)

            vertical=False
            if self.muAna.DSVcheck(detID): vertical=True 
            
            if vertical: 
                barposition=0.5*(A.x()+B.x()) 
                lam=(z-self.tw.pos.z())/self.tw.mom.z()
                tmp=ROOT.TVector3(self.tw.pos.x()+lam*self.tw.mom.x(), self.tw.pos.y()+lam*self.tw.mom.y(), self.tw.pos.z()+lam*self.tw.mom.z())             
                histname=f'pull_{key}'
                self.hists[histname].Fill(barposition - tmp.x(), nbars)  
            
            else: 
                barposition=0.5*(A.y()+B.y()) 
                lam=(z-self.tw.pos.z())/self.tw.mom.z()
                tmp=ROOT.TVector3(self.tw.pos.x()+lam*self.tw.mom.x(), self.tw.pos.y()+lam*self.tw.mom.y(), self.tw.pos.z()+lam*self.tw.mom.z())             
                histname=f'pull_{key}'
                self.hists[histname].Fill(barposition - tmp.y(), nbars)
         
    def FillAverageBarTime(self, hits, state=None):
        if not state: state=self.state
        
        for hit in hits:
            detID=hit.GetDetectorID()
            s,p,b=self.muAna.parseDetID(detID)
            if s==3: continue
            medians=self.muAna.GetMedianTime(hit, mode=state)
            if not medians: continue
            if len(medians)!=2: continue
            
            # zEx=self.zPos['MuFilter'][s*10+p]
            # lam=(zEx-self.tw.pos.z())/self.tw.mom.z()
            # Ex=ROOT.TVector3(self.tw.pos.x()+lam*self.tw.mom.x(), self.tw.pos.y()+lam*self.tw.mom.y(), self.tw.pos.z()+lam*self.tw.mom.z())
            self.MuFilter.GetPosition(detID, A, B)
            xpred=ROOT.TMath.Sqrt((A.x()-self.tw.Ex.x())**2 + (A.y()-self.tw.Ex.y())**2 + (A.z()-self.tw.Ex.z())**2)
            
            averagetime=sum(medians.values())/2 
            hist=self.hists[f'{s*10+p}_bar{b}_averagetime_{self.state}']
            hist.Fill(xpred, averagetime)
            
    def FillPlaneMultiplicity(self, hits):
        
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
                
    def FillTimingDiscriminantHists(self, hits):
        
        timingdiscriminant=self.muAna.GetTimingDiscriminant(hits, self.tw.MuFilter)
        if timingdiscriminant==-420:
            return
        
        self.hists['timingdiscriminant'].Fill(timingdiscriminant)
        if self.slopecut: self.hists['timingdiscriminant-slopecut'].Fill(timingdiscriminant)
        
        ds3haverage=self.muAna.GetDSHaverage(hits, mode='timingdiscriminant')
        if ds3haverage>-420: 
            us1averagetime = -1*(timingdiscriminant - ds3haverage)
            self.hists['averageDS3Htime'].Fill(ds3haverage)
            self.hists['averageUS1time'].Fill(us1averagetime)                   
                
    def FillFiredBarCorrelation(self, hits):
        
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
                    self.tw.MuFilter.GetPosition(detID, A, B)
                    docas[self.Getyresidual(detID)]=bar

                closest_bar=docas.pop(min(docas))
                for x in docas: self.hists[f'{10*subsystem+plane}_barcorrelation'].Fill(closest_bar, docas[x])

    def FilldeltaDSH(self, hits):
        res=self.muAna.GetDSHaverage(hits, mode='deltastations')
        if not res: return
        [self.hists[i].Fill(res[i]) for i in res]

        if self.tw.hasTrack: 
            self.deltads32vcuts['has track'].Fill(res['delta32'])
            self.deltads21vcuts['has track'].Fill(res['delta21'])
        if self.tw.passredchi2cut:
            self.deltads32vcuts['#chi^{2}_{#nu}'].Fill(res['delta32'])
            self.deltads21vcuts['#chi^{2}_{#nu}'].Fill(res['delta21'])
        if self.tw.passslopecut:
            self.deltads32vcuts['slope'].Fill(res['delta32'])
            self.deltads21vcuts['slope'].Fill(res['delta21'])

        # Fill times for DSH planes

    def FilldeltaSF(self, scifihits):
        res=self.muAna.GetScifiAverageTime(self.Scifi, scifihits, 'deltastations')
        [self.hists[i].Fill(res[i]) for i in res]

    def FillScifiDSresidual(self, hits):
        
        scifi_track=None
        for track in self.tw.M.Reco_MuonTracks:
            if track.GetUniqueID()==1:
                scifi_track=track 
                break 
        if scifi_track==None: return # No track in the Scifi

        scifi_fstate=scifi_track.getFittedState()
        scifi_fitStatus=scifi_track.getFitStatus()
        if not scifi_fitStatus.isFitConverged(): return
        scifi_pos = scifi_fstate.getPos()
        scifi_mom = scifi_fstate.getMom()
        
        for hit in hits:
            detID=hit.GetDetectorID()
            s,p,b=self.muAna.parseDetID(detID)
            if s!=2: continue
            
            zEx=self.zPos['MuFilter'][s*10+p]
            ds_lam=(zEx-self.tw.pos.z())/self.tw.mom.z()
            Ex = self.tw.Ex[p]
            scifi_lam=(zEx-scifi_pos.z())/scifi_mom.z()
            scifi_Ex=ROOT.TVector3(scifi_pos.x()+scifi_lam*scifi_mom.x(), scifi_pos.y()+scifi_lam*scifi_mom.y(), scifi_pos.z()+scifi_lam*scifi_mom.z())

            self.hists[f'DSxvScifix_plane{p}_bar{b}'].Fill(Ex.x() - scifi_Ex.x())
            
    def FillChannelRateHists(self):
        hits=self.M.eventTree.Digi_MuFilterHits
        for hit in hits:
            detID=hit.GetDetectorID()
            s,p,b=self.muAna.parseDetID(detID)

            if s==3: continue
            for x in hit.GetAllTimes():
                SiPM, clock = x
                channel=self.muAna.GetSiPMNumberInSystem_PCBbyPCB(detID, SiPM)
                self.hists[f'channelhitrate-{self.subsystemdict[s]}'].Fill(channel)

                if self.tw.hasTrack:

                    td=self.muAna.GetTimingDiscriminant()
                    if not self.TimingDiscriminantCut(td):
                        self.hists[f'channelhitrate-{self.subsystemdict[s]}-tdcut'].Fill(channel)

                    if not self.yresidual3(detID): 
                        self.hists[f'channelhitrate-{self.subsystemdict[s]}-dycut'].Fill(channel)

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

    def Filltds0vcuts(self):

        if self.tw.hasTrack: 
            self.tds0vcuts['has track'].Fill(self.tw.reft)
        if self.tw.passredchi2cut:
            self.tds0vcuts['#chi^{2}_{#nu}'].Fill(self.tw.reft)
        if self.tw.passslopecut:
            self.tds0vcuts['slope'].Fill(self.tw.reft)

    def WriteOutHistograms(self):

        outpath=f'{self.outpath}splitfiles/run{self.runNr}/selectioncriteria/'
        path_obj=Path(outpath)
        path_obj.mkdir(parents=True, exist_ok=True)
        outfile=f'selectioncriteria_{self.options.nStart}.root'
        f=ROOT.TFile.Open(outpath+outfile, 'recreate')
        
        additionalkeys=['DSxvScifix']

        f.WriteObject(self.tw.hists[f'reft-{self.tw.refsysname}'], f'reft-{self.tw.refsysname}', 'kOverwrite')

        for histname in self.hists:
            if not any( [histname.find(additionalkey)>-1 for additionalkey in additionalkeys] ):
                
                hist=self.hists[histname]
                f.WriteObject(hist, hist.GetName(), 'kOverwrite')
            else:
                for additionalkey in additionalkeys:
                    if histname.find(additionalkey)>-1:
                        if not hasattr(f, additionalkey): residualsfolder=f.mkdir(additionalkey)
                        else: residualsfolder=f.Get(additionalkey)
                        residualsfolder.cd()
                        hist=self.hists[histname]
                        hist.Write(histname, 2) # The 2 means it will overwrite a hist of the same name
        
        # Make 2d hist comparing cut impact on tds0
        self.hists['tds0vcuts']=ROOT.TH2D('tds0vcuts', 'Effect of different cuts on the t_{0}^{DS} distribution;Cut name;t^{DS}_{0} [ns]', len(self.cuteffect), 0, len(self.cuteffect), 100, 0, 25)
        self.hists['deltads32vcuts']=ROOT.TH2D('deltads32vcuts', 'Effect of different cuts on #Delta (t_{DS3}, t_{DS2}) distribution;Cut name;#Delta (t_{DS3}, t_{DS2}) [ns]', len(self.cuteffect), 0, len(self.cuteffect), 80, -10, 10)
        self.hists['deltads21vcuts']=ROOT.TH2D('deltads21vcuts', 'Effect of different cuts on #Delta (t_{DS2}, t_{DS1}) distribution;Cut name;#Delta (t_{DS2}, t_{DS1}) [ns]', len(self.cuteffect), 0, len(self.cuteffect), 100, 0, 25)        

        for i, (key, value) in enumerate(self.cuteffect.items(), start=1):  # Bin numbers start at 1
            self.hists['cutflow'].GetXaxis().SetBinLabel(i, key)  # Set bin label
            self.hists['tds0vcuts'].GetXaxis().SetBinLabel(i, key)  # Set bin label
            self.hists['deltads32vcuts'].GetXaxis().SetBinLabel(i, key)  # Set bin label
            self.hists['deltads21vcuts'].GetXaxis().SetBinLabel(i, key)  # Set bin label
            self.hists['cutflow'].SetBinContent(i, value)  # Set bin content

            cutdist = self.tds0vcuts[key]
            for xbin in range(cutdist.GetNbinsX()):
                bincontent = cutdist.GetBinContent(xbin)
                self.hists['tds0vcuts'].SetBinContent(i, xbin, bincontent)

            cutdistd32 = self.deltads32vcuts[key]
            cutdistd21 = self.deltads21vcuts[key]
            for xbin in range(cutdistd32.GetNbinsX()):
                d32_bincontent = cutdistd32.GetBinContent(xbin)
                d21_bincontent = cutdistd21.GetBinContent(xbin)
                self.hists['deltads32vcuts'].SetBinContent(i, xbin, d32_bincontent)
                self.hists['deltads21vcuts'].SetBinContent(i, xbin, d21_bincontent)
                
        f.WriteObject(self.hists['cutflow'], 'cutflow', 'kOverwrite')
        f.WriteObject(self.hists['tds0vcuts'], 'tds0vcuts', 'kOverwrite')

        f.WriteObject(self.hists['deltads32vcuts'], 'deltads32vcuts', 'kOverwrite')
        f.WriteObject(self.hists['deltads21vcuts'], 'deltads21vcuts', 'kOverwrite')
        
        f.Close()
        print(f'{len(self.hists)} histograms saved to {outpath}{outfile}')
        
    def WriteOutToTestFile(self):       
        outfile=f'testing_{self.options.runNumber}_{self.options.nStart}_{self.options.nEvents}.root'
        f=ROOT.TFile.Open(outfile, 'update')
        
        # mainhistkeys=['dy', 'slopes', 'nSiPMs', 'timingdiscriminant', 'trackfitting', 'reft', ]
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