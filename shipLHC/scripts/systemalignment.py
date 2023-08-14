import ROOT, os, csv, subprocess, cProfile, io, pstats
from argparse import ArgumentParser

from AnalysisFunctions import Analysis
from pathlib import Path
from array import array

A, B, locA, locB = ROOT.TVector3(), ROOT.TVector3(), ROOT.TVector3(), ROOT.TVector3()

# class SystemAlignment(ROOT.FairTask):
class SystemAlignment(object):

    def __init__(self, options, tw):
       
        # self.state='corrected'
        # options.state='corrected'
        self.options=options
        self.tw=tw
        self.runNr = tw.runNr 
        self.muAna = tw.muAna
        self.timealignment=tw.timealignment

        ### If no time-walk correction run is provided. Set the default correction run depending on time alignment of the data set


        self.afswork=tw.afswork
        self.outpath=tw.outpath

        self.subsystemdict={1:'Veto', 2:'US', 3:'DS'}
        self.nchs={1:224, 2:800}

        self.systemAndPlanes = {1:2,2:5,3:7}
        self.systemAndBars={1:7,2:10,3:60}
        self.systemAndChannels={1:[8,0],2:[6,2],3:[1,0]}
        self.sdict={0:'Scifi',1:'Veto',2:'US',3:'DS'}
        self.zPos=tw.zPos
        self.cutdists=self.muAna.GetCutDistributions(self.runNr, ('dy', 'timingdiscriminant'))

        self.freq=160.316E6
        self.TDC2ns=1E9/self.freq

        self.xrefs=tw.xrefs
        self.sides=('left', 'right')

        self.hists=tw.hists

        # self.gauss=ROOT.TF1('mygauss', 'abs([0])*abs([4])/(abs([2])*sqrt(2*pi))*exp(-0.5*((x-[1])/[2])**2)+abs([3])', 4)
        # for idx,param in enumerate(('Constant', 'Mean', 'Sigma', 'y-offset', 'bin-width')): self.gauss.SetParName(idx, param)

        self.sigmatds0=0.263 # ns 

    def FillSiPMHists(self, hit):

        detID=hit.GetDetectorID()

        # qdcs, clocks = hit.GetAllSignals(), hit.GetAllTimes()

        alignedtimes=self.muAna.GetCorrectedTimes(hit, self.tw.pred, mode='unaligned')

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
            if self.options.debug: print(f'{fixed_ch}, {self.tw.TDS0 - time - alignmentparameter[0] }')

            self.hists[SiPMtime].Fill(self.tw.TDS0 - time - alignmentparameter[0])

    def FillBarTimeHists(self, hit):

        """
        I will keep the Analysis class methods to return the 
        SiPM times relative to the event-t0. In the methods in this SystemAlignment class
        I will shift them to be relative to tds0.
        """
        
        ### medians returns a dictionary: {'left': float, 'right':float}
        # medians=self.muAna.GetMedianTime(hit, mode='alignment')
        alignedtimes=self.muAna.GetAlignedTimes(hit, self.tw.pred)
        
        if medians==-999: return
        if not medians: 
            return 

        averagetime=1/2 * sum([self.tw.TDS0-i for i in medians.values()])

        deltatime=self.tw.TDS0-medians['left'] - self.tw.TDS0-medians['right']
        
        detID=hit.GetDetectorID()
        subsystem, plane, bar= self.muAna.parseDetID(detID)
        
        averagebartimehistname=f'{10*subsystem+plane}_bar{bar}_averagetime'
        if not averagebartimehistname in self.hists:
            title=self.subsystemdict[subsystem]+' plane '+str(plane+1)+' bar '+str(bar+1)+' average of median '+self.state+' times from each side as a function of x_{predicted};x_{predicted} [cm];'+'#frac{1}{2}#times(t^{tw corr}_{left}+t^{tw corr}_{right}) [ns];Counts'
            self.hists[averagebartimehistname]=ROOT.TH2F(averagebartimehistname, title, 100, 0, 100, 200, -10, 10)
        self.hists[averagebartimehistname].Fill(self.pred, self.tw.TDS0-averagetime)
        
        deltabartimehistname=f'{10*subsystem+plane}_bar{bar}_deltatime'
        if not deltabartimehistname in self.hists:
            title=self.subsystemdict[subsystem]+' plane '+str(plane+1)+' bar '+str(bar+1)+' difference of median '+self.state+' times from each side as a function of x_{predicted};x_{predicted} [cm];'+'t^{tw corr}_{left} - t^{tw corr}_{right} [ns];Counts'
            self.hists[deltabartimehistname]=ROOT.TH2F(deltabartimehistname, title, 100, 0, 100, 200, -10, 10)
        self.hists[deltabartimehistname].Fill(self.tw.pred, deltatime)

        for side in ('left', 'right'):
            averagebarsidetimehistname=f'{10*subsystem+plane}_bar{bar}_{side}sidetime'
            if not averagebarsidetimehistname in self.hists:
                title=self.subsystemdict[subsystem]+' plane '+str(plane+1)+' bar '+str(bar+1)+' median '+self.state+' times from '+side+' side as a function of x_{predicted};x_{predicted} [cm];'+'Median of t^{tw corr}_{'+side+'} [ns];Counts'
                self.hists[averagebarsidetimehistname]=ROOT.TH2F(averagebarsidetimehistname, title, 100, 0, 100, 200, -10, 10)
            self.hists[averagebarsidetimehistname].Fill(self.tw.pred, self.tw.TDS0- medians[side])

    def XTHists(self, hit):

        clocks, qdcs =hit.GetAllTimes(), hit.GetAllSignals()
        times={'left':{}, 'right':{}}
        detID=hit.GetDetectorID()
        s,p,b = self.muAna.parseDetID(detID)

        for i in clocks: 
            SiPM, clock=i
            fixed_ch=self.muAna.MakeFixedCh((s,p,b,SiPM))
            qdc=self.muAna.GetChannelVal(SiPM, qdcs)
            correctedtime=self.muAna.MuFilterCorrectedTime(fixed_ch, qdc, clock, self.tw.pred)
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
                    ti, tj= self.tw.TDS0-times[side][i],self.tw.TDS0-times[side][j]
                    self.hists[xthistname].Fill(ti, tj)

    def WriteOutHistograms(self):

        outfilename=f'{self.outpath}splitfiles/run{self.runNr}/SystemAlignment/SystemAlignment_{self.options.nStart}.root'
        if os.path.exists(outfilename): outfile=ROOT.TFile.Open(outfilename, 'recreate')
        else: outfile=ROOT.TFile.Open(outfilename, 'create')

        additionalkeys=['averagetime', 'deltatime', 'sidetime', 'AlignedSiPMtime', 'timingxt']
        for h in self.hists:
            if h=='TDS0':
                outfile.WriteObject(self.hists[h], self.hists[h].GetName(),'kOverwrite')

            if len(h.split('_'))==3:
                for additionalkey in additionalkeys:
                    if h.find(additionalkey)==-1: continue
                    planekey, bar, key = h.split('_')
                    
                    hist=self.hists[h]
                    if not hasattr(outfile, additionalkey): folder=outfile.mkdir(additionalkey)
                    else: folder=outfile.Get(additionalkey)
                    folder.cd()
                    hist.Write(h, 2) # The 2 means it will overwrite a hist of the same name                

        outfile.Close()
        print(f'{len(self.hists)} histograms saved to {outfilename}')        

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
                      