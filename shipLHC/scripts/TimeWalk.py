import ROOT, os, csv, subprocess, cProfile, io, pstats
from argparse import ArgumentParser
import rootUtils as ut
import AnalysisFunctions as muAna
import SndlhcGeo, SndlhcTracking
from pathlib import Path

deciles=[i/10 for i in range(11)]
v1, v2 = ROOT.TVector3(), ROOT.TVector3()

class TimeWalk(ROOT.FairTask):

    #def Init(self, options, runtw):
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
        self.outpath=self.afswork
        self.EventNumber=-1

        self.systemAndPlanes = {1:2,2:5,3:7}
        self.systemAndBars={1:7,2:10,3:60}
        self.systemAndChannels={1:[8,0],2:[6,2],3:[1,0]}
        self.sdict={0:'Scifi',1:'Veto',2:'US',3:'DS'}
        self.zPos=self.M.zPos
        self.cutdists=muAna.GetCutDistributions(self.runNr, 'dy', options.nStations)
        iterationdict={'zeroth':0, 'ToF':1, 'TW':2}
        self.iteration=iterationdict[options.mode]

        self.freq=160.316E6
        self.TDC2ns=1E9/self.freq

        self.largeSiPMmap={0:0 ,1:1 ,3:2 ,4:3 ,6:4 ,7:5}
        self.verticalBarDict={0:1, 1:3, 2:5, 3:6}
        verticalPlanes=list(self.verticalBarDict.values())
        self.xref=42. # To be set as USbarlength/2.
        self.correctionfunction=lambda ps, qdc, limit: 1/sum( [ ps[i]*qdc**i for i in range(len(ps)) ] ) 
        self.hists=self.M.h

        self.hists['nStations'] = ROOT.TH1F('nStations', 'Number of stations used in track fit;Number of stations used;Counts', 10, 0, 10)
        self.hists['nStations'].Fill(options.nStations)        
        
        if options.mode=='TW':
            if not options.CorrectionType: self.CorrectionType=4

            self.correctionparams=lambda ps : [y for x,y in enumerate(ps) if x%2==0]
            self.CorrectionType=options.CorrectionType
            if options.CorrectionType==1: self.correctionfunction = lambda ps, qdc : 1/sum( [ ps[i]*qdc**i for i in range(len(ps)) ] ) 
            elif options.CorrectionType==5: self.correctionfunction = lambda ps, qdc : ps[0]*ROOT.Log(ps[1]*qdc)+ps[2]
            elif options.CorrectionType==4: self.correctionfunction = lambda ps, qdc : ps[4]+ ps[3](qdc-ps[0])/( ps[0] + ps[1]*(qdc-ps[0]) + ps[2]*(qdc-ps[0])*(qdc-ps[0]) ) 

    def GetEntries(self):
        return self.eventTree.GetEntries()

    def ExecuteEvent(self, event):

        hists=self.hists
        
        Reco_MuonTracks=self.M.Reco_MuonTracks
        theTrack=Reco_MuonTracks[0]
        fstate=theTrack.getFittedState()
        pos=fstate.getPos()
        mom=fstate.getMom()

        # Returns bool for if the track object is in the acceptance for Veto and US respectively.
        InVeto, InUS=(muAna.InAcceptance(pos, mom, i, self.MuFilter, self.zPos) for i in (1,2))
        if not InUS: return 0

        # Use the correct subsystems for 1 bar / plane cut
        if InVeto: tmp=(1,2)
        else: tmp=(2,)
        if not muAna.OneHitPerSystem(event.Digi_MuFilterHits, tmp): return 0

        # Get DS event t0
        DST0cc=muAna.GetDSH_average(event.Digi_MuFilterHits)
        if not 'DST0' in self.hists:
           self.hists['DST0']=ROOT.TH1F('DST0','DSH average', 100, 0, 25)
        self.hists['DST0'].Fill(DST0cc*6.25)
        DST0=self.DST0cut(DST0cc)
        if DST0==0: return 0

        ### Slope cut
        slopeX, slopeY = mom.x()/mom.z(), mom.y()/mom.z()
        if not TimeWalk.slopecut(mom): return 0

        for hit in event.Digi_MuFilterHits:
            detID=hit.GetDetectorID()
            s,p,b=muAna.parseDetID(detID)
            if s==3: continue
            if not self.yresidual(detID,pos,mom): continue

            zEx=self.zPos['MuFilter'][s*10+p]
            lam=(zEx-pos.z())/mom.z()
            Ex=ROOT.TVector3(pos.x()+lam*mom.x(), pos.y()+lam*mom.y(), pos.z()+lam*mom.z())

            xpred=self.xpred(v1,Ex)
            xpred_R=-xpred

            channels_t=hit.GetAllTimes()
            channels_qdc=hit.GetAllSignals()
            for channel in channels_t:
                SiPM,clock=channel
                qdc=muAna.GetChannelVal(SiPM, channels_qdc)
                if qdc==-999.: continue
                fixed=(s,p,b,SiPM)
                fixed_ch=muAna.MakeFixedCh(fixed)

                if self.options.mode=='zeroth': self.zeroth(fixed_ch,xpred,clock,qdc,DST0)
                elif self.options.mode=='ToF': self.ToF(fixed_ch,xpred, clock, qdc, DST0)
                elif self.options.mode=='TW': self.TW(fixed_ch, xpred, clock, qdc, DST0)

    def zeroth(self,fixed_ch,xpred,clock,qdc,DST0):
        hists=self.hists
        
        detID, SiPM=( int(fixed_ch.split('_')[i]) for i in range(2) ) 

        correctedtime=clock*6.25
        t_rel=DST0-correctedtime
        dtvxpred=f'dtvxpred_{fixed_ch}_iteration{self.iteration}'
        if not dtvxpred in hists:
           title='Uncorrected T_{0}^{DS}-t_{'+str(SiPM)+'} v x-position;x_{predicted} [cm];T_{0}^{DS}-t_{'+str(SiPM)+'} [ns]'
           hists[dtvxpred]=ROOT.TH2F(dtvxpred,title,110,-10,100, 160, -20, 20.)

        SiPMtime=f'tSiPM_{fixed_ch}_iteration{self.iteration}'
        if not SiPMtime in hists:
           hists[SiPMtime]=ROOT.TH1F(SiPMtime,f'Measured SiPM time {fixed_ch};SiPM time [ns];Counts', 100, 0, 25)
        attlen=f'attlen_{fixed_ch}_iteration{self.iteration}'
        if not attlen in hists:
            title='Predicted position against QDC_{'+str(SiPM)+'} '+fixed_ch+';x_{predicted} [cm];QDC_{'+str(SiPM)+'} [a.u]'
            hists[attlen]=ROOT.TH2F(attlen,title, 110, -10, 110, 200, 0., 200)
        # trackposveto0=f'DS muon track position at veto 0;Track x [cm];Track y [cm],'    
        
        self.hists[dtvxpred].Fill(xpred,t_rel)
        self.hists[SiPMtime].Fill(correctedtime)
        self.hists[attlen].Fill(xpred, qdc)

    def ToF(self, fixed_ch, xpred, clock, qdc, DST0):
        hists=self.hists

        cdata=muAna.Getcscint(self.runNr, fixed_ch, self.iteration-1)
        if cdata==-999.: return 0
        SiPM=int(fixed_ch.split('_')[-1])
        correctedtime=muAna.correct_ToF(SiPM, clock, xpred, cdata, self.xref)[1]
        name=f'dtvqdc_{fixed_ch}_iteration{self.iteration}'
        if not name in hists:
            title='Uncorrected T_{0}^{DS}-t_{'+str(SiPM)+'} v QDC_{'+str(SiPM)+'};QDC_{'+str(SiPM)+'} [a.u];T_{0}^{DS}-t_{'+str(SiPM)+'} [ns]'
            hists[name]=ROOT.TH2F(name, title, 200, 0, 200, 160, -20, 20)
        t_rel=DST0-correctedtime
        self.hists[name].Fill(qdc,t_rel)

    def TW(self, fixed_ch, xpred, clock, qdc, DST0):
        hists=self.hists
        
        SiPM=int(fixed_ch.split('_')[-1])
        time=clock*6.25
        paramsAndErrors,correctionlimits = muAna.GetPolyParams(self.runNr, fixed_ch, self.CorrectionType, 0)
        chi2pNDF = muAna.Getchi2pNDF(self.runNr, fixed_ch, 0)
        cdata = muAna.Getcscint(self.runNr, fixed_ch, 0)
        
        if paramsAndErrors==-999. or cdata==-999.: return 0
        polyparams = self.correctionparams(paramsAndErrors)

        if chi2pNDF>2:return 0
       
        # Only if qdc is less than the high QDC limit from the fit: apply the correction
        if qdc <= correctionlimit[1]:
            if len(polyparams)==4:
                polycorrection = self.correctionfunction(polyparams[0:-1], qdc)
                correctedtime = time + polycorrection+polyparams[-1]
            elif len(polyparams)==3:
                polycorrection = self.correctionfunction(polyparams, qdc)
                correctedtime = time + polycorrection
        
            correctedtime2 = muAna.correct_ToF(SiPM, clock, xpred, cdata, self.xref)[1] + polycorrection
        
        else: 
            correctiontime=time
            correctedtime2=muAna.correct_ToF(SiPM, clock, xpred, cdata, self.xref)[1]
        
        t_rel=DST0 - correctedtime
        t2_rel=DST0 - correctedtime2
       
        if self.options.debug:
            ### Debugging TW correction application
            print(f'-----{fixed_ch}-----\ntime: {time}\npolycorrection: {polycorrection}')
            if len(polyparams)==4: print(f'\noffset:{polyparams[-1]}')
            print('---------------\n')

        name=f'dtvxpred_{fixed_ch}_iteration{self.iteration}'
        if not name in hists:
           title='Corrected T_{0}^{DS}-t_{'+str(SiPM)+'} v x-position;x_{predicted} [cm];T_{0}^{DS}-t_{'+str(SiPM)+'} [ns]'
           hists[name]=ROOT.TH2F(name,title,110,-10,100, 160, -20, 20.)
        
        name2=f'dtvqdc_{fixed_ch}_iteration{self.iteration}'
        if not name2 in hists:
            title='Corrected T_{0}^{DS}-t_{'+str(SiPM)+'} v QDC_{'+str(SiPM)+'};QDC_{'+str(SiPM)+'} [a.u];T_{0}^{DS}-t_{'+str(SiPM)+'} [ns]'
            hists[name2]=ROOT.TH2F(name2, title, 200, 0, 200, 160, -20, 20)

        hists[name].Fill(xpred,t_rel)
        hists[name2].Fill(qdc, t2_rel)

    def FillDST0(self):
        pass 


    def WriteOutHistograms(self):
        
        for h in self.M.h:
            if len(h.split('_'))==4:
                histkey,detID,SiPM,iterationNumber=h.split('_')
                fixed_ch='_'.join((detID,SiPM))
                hist=self.M.h[h]
                outpath=f'{self.afswork}/splitfiles/run{self.runNr}/{fixed_ch}/'
                path_obj=Path(outpath)
                path_obj.mkdir(parents=True, exist_ok=True)
                
                outfile=f'timewalk_{fixed_ch}_{self.options.nStart}.root'
                f=ROOT.TFile.Open(outpath+outfile, 'update')
                f.WriteObject(hist, hist.GetName(), 'kOverwrite')
                if self.options.mode=='zeroth':
                    f.WriteObject(self.hists['DST0'], 'DST0', 'kOverwrite')
                f.Close()
        print(f'{len(self.M.h)} histograms saved to {self.afswork}/run{self.runNr}')

    def yresidual(self,detID,pos,mom): 
        s,p,b=muAna.parseDetID(detID) 
        key=10*s+p 
        z=self.zPos['MuFilter'][key] 
        lam=(z-pos.z())/mom.z() 
        Ex=ROOT.TVector3(pos.x()+lam*mom.x(), pos.y()+lam*mom.y(), pos.z()+lam*mom.z()) 
        self.MuFilter.GetPosition(detID,v1,v2) 
        dy=Ex.y()-v1.y() # Needs to be checked 
        dy_min, dy_max = TimeWalk.dycut(self.cutdists[f'dy_{key}_{self.options.nStations}stations']) 
        #if dymin or dy 
        if dy>dy_max or dy<dy_min: return 0 
        else: return 1 
    
    @staticmethod
    def dycut(hist, nsig=1):    
        dymin=hist.GetMean()-nsig*hist.GetStdDev()
        dymax=hist.GetMean()+nsig*hist.GetStdDev()
        return dymin, dymax

    @staticmethod
    def slopecut(mom, slopecut=0.4):
        slopeX, slopeY=mom.x()/mom.z(), mom.y()/mom.z()
        if abs(slopeX)>slopecut or abs(slopeY)>slopecut: return 0
        else: return 1

    def DST0cut(self,DST0):
        if DST0==-999. or DST0==-998. or DST0==-6237.5: return 0
        DST0ns=DST0*6.25
        if DST0ns<13.75 or DST0ns>15.45: return 0
        return DST0ns 
    
    def xpred(self,v1,Ex):
        return ROOT.TMath.Sqrt((v1.x()-Ex.x())**2+(v1.y()-Ex.y())**2+(v1.z()-Ex.z())**2)


