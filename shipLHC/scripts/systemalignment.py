#!/usr/bin/env python
import ROOT, os

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

        if options.XT:
            from itertools import combinations
            self.vetocombinations={'left':list(combinations(range(8), 2)), 'right':list(combinations(range(8, 16), 2))}
            self.UScombinations={'left':list(combinations([0,1,3,4,6,7], 2)), 'right':list(combinations([8,9,11,12,14,15], 2))}

            self.SiPMcombinations = {1:{'left':list(combinations(range(8), 2)), 'right':list(combinations(range(8, 16), 2))}, 
            2:{'left':list(combinations([0,1,3,4,6,7], 2)), 'right':list(combinations([8,9,11,12,14,15], 2))}
            }
        
        self.sigmatds0=0.263 # ns 

    def FillSiPMHists(self, hit):
        detID=hit.GetDetectorID()

        correctedtimes=self.muAna.GetCorrectedTimes(hit, self.tw.pred, mode='aligned')

        for ch in correctedtimes:
            SiPM, time=ch
            fixed_ch=f'{detID}_{SiPM}'

            ReadableDetID=self.muAna.MakeHumanReadableFixedCh(fixed_ch)
            SiPMtime=f'AlignedSiPMtime_{fixed_ch}'
            if not SiPMtime in self.hists:
                title=f'DS horizontal average relative, TW + ToF corrected + DS aligned SiPM time'
                splittitle='#splitline{'+ReadableDetID+'}{'+title+'}'
                axestitles='t^{DS}_{0} - t_{SiPM '+str(SiPM)+'}^{tw corr} - d_{SiPM '+str(SiPM)+'} [ns];Counts'
                fulltitle=splittitle+';'+axestitles
                # if self.timealignment=='old': self.hists[SiPMtime]=ROOT.TH1F(SiPMtime,fulltitle, 400, -10, 10)
                self.hists[SiPMtime]=ROOT.TH1F(SiPMtime,fulltitle, 400, -10, 10)
                # else: self.hists[SiPMtime]=ROOT.TH1F(SiPMtime,fulltitle, 2000, -5, 5)

            d=self.muAna.alignmentparameters[fixed_ch]

            self.hists[SiPMtime].Fill(time)

    def FillBarHists(self, hit):
        # if not self.tw.passslopecut: return
        
        detID=hit.GetDetectorID()
        s, p, b= self.muAna.parseDetID(detID)

        # correctedtimes=self.muAna.GetCorrectedTimes(hit, self.tw.pred, mode='unaligned')
        correctedtimes = self.muAna.GetCorrectedTimes(hit, x=0, mode='aligned')
        qdcs = hit.GetAllSignals()
        aligned_times={'left':{}, 'right':{}}
        times={'left':{}, 'right':{}} # Only apply time walk correction

        for i in correctedtimes: 
            SiPM, correctedtime = i
            fixed_ch=self.muAna.MakeFixedCh((s,p,b,SiPM))

            if not fixed_ch in self.muAna.alignmentparameters: continue
            d=self.muAna.alignmentparameters[fixed_ch]

            qdc=self.muAna.GetChannelVal(SiPM, qdcs)
            side=self.muAna.GetSide(f'{detID}_{SiPM}')
            # aligned_times[side][SiPM] = self.tw.reft - correctedtime - d[0]
            aligned_times[side][SiPM] = correctedtime
            # times[side][SiPM] = correctedtime - d[0]

        if any([len(aligned_times[i])==0 for i in ('left', 'right')]): return
        averages={side:self.tw.reft - (sum(aligned_times[side].values()) / len(aligned_times[side])) for side in aligned_times}
        fastest = {side:self.tw.reft - min(aligned_times[side].values()) for side in aligned_times}
        averagetime = 1/2 * sum(averages.values())

        baraveragecscint_left, baraveragecscint_right = self.muAna.GetBarAveragecscint(self.runNr, detID, 'corrected')
        baraveragecscint = 0.5* (baraveragecscint_left[0] + baraveragecscint_right[0])

        delta_averagetime = averages['right'] - averages['left']
        delta_fastesttime = fastest['right'] - fastest['left']

        # Using passing muons, bar length = 2*x_L / c dt_{L,R}
        # where x_L is the distance from the left edge of the bar to the extrapolated track
        
        # self.pred is set for each hit in the TimeWalk task
        x='yAgreement' if self.tw.trackrelated else 'notyAgreement'
        xL = f'xL-{x}_plane{p}'
        if not xL in self.hists:
            title='x_{L};x_{L} [cm];Counts'
            self.hists[xL]=ROOT.TH1F(xL,title,80, 0, 80)
        self.hists[xL].Fill(self.tw.pred)

        # Looking at differences depending on whether the average or fastest times are used
        for timer in ('fastest', 'average'):
            dt = f'dt-{x}_plane{p}_{timer}'
            if not dt in self.hists:
                title='#Delta t(left, right);#Delta t_{L, R} [ns];Counts'
                self.hists[dt]=ROOT.TH1F(dt,title,50, -5, 5)
            
            barycentrex_hist = f'xbarycentre-{x}_plane{p}_{timer}'
            if not barycentrex_hist in self.hists:
                title='Barycentre in x, #frac{c_{scint}}{2}(t_{right}+t_{left});#frac{c_{scint}}{2}(t_{right}+t_{left}) [cm];Counts'
                self.hists[barycentrex_hist] = ROOT.TH1F(barycentrex_hist, title, 160, -90, 60)
               
            sumt = f'sumt-{x}_plane{p}_{timer}'
            if not sumt in self.hists:
                title = 'corrected and aligned t_{right} + t_{left};t_{right} + t_{left} [ns];Counts'

            lamx_hist = f'lambda-{x}_plane{p}_{timer}'
            if not lamx_hist in self.hists:
                title='Deposition width c_{scint}(t_{right} + t_{left});c_{scint}(t_{right} + t_{left}) [cm];Counts'
                self.hists[lamx_hist] = ROOT.TH1F(lamx_hist, title, 160, -90, 60)

            barycentrex_hist = f'xbarycentre-{x}_plane{p}_{timer}'
            if not barycentrex_hist in self.hists:
                title='Barycentre in x, #frac{c_{scint}}{2}(t_{right}+t_{left});#frac{c_{scint}}{2}(t_{right}+t_{left}) [cm];Counts'
                self.hists[barycentrex_hist] = ROOT.TH1F(barycentrex_hist, title, 160, -90, 60)
                

            if timer=='fastest':        
                self.hists[dt].Fill(delta_fastesttime)
                lamx = baraveragecscint*delta_fastesttime
                barycentrex = 0.5*baraveragecscint*(fastest['right'] + fastest['left'])
                self.hists[lamx_hist].Fill(lamx)
                self.hists[barycentrex_hist].Fill(barycentrex)
            
            elif timer=='average':
                self.hists[dt].Fill(delta_averagetime)
                lamx = baraveragecscint*delta_averagetime
                barycentrex = 0.5*baraveragecscint*(averages['right'] + averages['left'])
                self.hists[lamx_hist].Fill(lamx)
                self.hists[barycentrex_hist].Fill(barycentrex)
        
        averagebartimehistname=f'averagetime_{detID}_aligned'
        if not averagebartimehistname in self.hists:
            title=self.subsystemdict[s]+' plane '+str(p+1)+' bar '+str(b+1)+' average of aligned times from each side as a function of x_{predicted};x_{predicted} [cm];'+'#frac{1}{2}#times(t^{tw corr}_{left}+t^{tw corr}_{right}) [ns];Counts'
            self.hists[averagebartimehistname]=ROOT.TH2F(averagebartimehistname, title, 100, 0, 100, 2000, -5, 5)
        self.hists[averagebartimehistname].Fill(self.tw.pred, averagetime)
        
        deltabartimehistname=f'deltatime_{detID}_aligned'
        if not deltabartimehistname in self.hists:
            title=self.subsystemdict[s]+' plane '+str(p+1)+' bar '+str(b+1)+' difference of average aligned time from each side as a function of x_{predicted};x_{predicted} [cm];'+'t^{tw corr}_{left} - t^{tw corr}_{right} [ns];Counts'
            self.hists[deltabartimehistname]=ROOT.TH2F(deltabartimehistname, title, 100, 0, 100, 2000, -5, 5)
        self.hists[deltabartimehistname].Fill(self.tw.pred, delta_fastesttime)

        for side in ('left', 'right'):
            averagebarsidetimehistname=f'sidetime_{detID}-{side}_aligned'
            if not averagebarsidetimehistname in self.hists:
                title=self.subsystemdict[s]+' plane '+str(p+1)+' bar '+str(b+1)+' average aligned time from '+side+' side as a function of x_{predicted};x_{predicted} [cm];'+ side+' side average of t^{tw corr}_{'+side+'} [ns];Counts'
                self.hists[averagebarsidetimehistname]=ROOT.TH2F(averagebarsidetimehistname, title, 100, 0, 100, 2000, -5, 5)
            self.hists[averagebarsidetimehistname].Fill(self.tw.pred, averages[side])

    def XTHists(self, hit):

        correctedtimes, qdcs = self.muAna.GetCorrectedTimes(hit, self.tw.pred, mode='aligned'), hit.GetAllSignals()
        times={'left':{}, 'right':{}}
        detID=hit.GetDetectorID()
        s,p,b = self.muAna.parseDetID(detID)

        for i in correctedtimes: 
            SiPM, correctedtime = i
            fixed_ch=self.muAna.MakeFixedCh((s,p,b,SiPM))

            if not fixed_ch in self.muAna.alignmentparameters: continue
            d=self.muAna.alignmentparameters[fixed_ch]

            side=self.muAna.GetSide(f'{detID}_{SiPM}')
            times[side][SiPM] = correctedtime

        xrefthistname=f'timingxt_reft_{self.subsystemdict[s]}'
        if not xrefthistname in self.hists:
            title='Timing correlation between all '+self.subsystemdict[s]+' SiPMs and t_{0}^{DS}'
            axestitles='t_{0}^{DS} - t_{SiPM}^{tw corr} - d_{SiPM} [ns];t_{0}^{DS} [ns];Counts'
            fulltitle=title+';'+axestitles
            if self.timealignment=='old': self.hists[xrefthistname]=ROOT.TH2F(xrefthistname,fulltitle, 150, -5, 20, 125, 0, 25)
            else: self.hists[xrefthistname]=ROOT.TH2F(xrefthistname,fulltitle, 150, -5, 20, 125, 0, 25)

        for side in times:
            for combination in self.SiPMcombinations[s][side]:
                i, j=combination
                if not (i in times[side] and j in times[side]): continue
                
                time_i, time_j = times[side][i], times[side][j]
                xthistname=f'timingxt_{detID}_SiPMs{i}-{j}'
                if not xthistname in self.hists:
                    title=f'Timing correlation, {self.subsystemdict[s]} plane {p+1} bar {b+1}'
                    subtitle=f'{side} SiPMs {i+1} and {j+1}'
                    splittitle='#splitline{'+title+'}{'+subtitle+'}'
                    axestitles='t_{0}^{DS} - t_{SiPM '+str(i)+'}^{tw corr} [ns];t_{0}^{DS} - t_{SiPM '+str(j)+'}^{tw corr} [ns];Counts'
                    fulltitle=splittitle+';'+axestitles
                    if self.timealignment=='old': self.hists[xthistname]=ROOT.TH2F(xthistname,fulltitle, 150, -5, 20, 150, -5, 20)
                    else: self.hists[xthistname]=ROOT.TH2F(xthistname,fulltitle, 150, -5, 20, 150, -5, 20)

                self.hists[xthistname].Fill(time_i, time_j)
                self.hists[xrefthistname].Fill(time_i, self.tw.reft)

    def WriteOutHistograms(self):

        outfilename=f'{self.outpath}splitfiles/run{self.runNr}/SystemAlignment/SystemAlignment_{self.options.nStart}.root'
        if os.path.exists(outfilename): outfile=ROOT.TFile.Open(outfilename, 'recreate')
        else: outfile=ROOT.TFile.Open(outfilename, 'create')

        additionalkeys=['averagetime', 'deltatime', 'sidetime', 'AlignedSiPMtime', 'timingxt',
                        'xL', 'dt', 'sumt', 'lambda', 'xbarycentre']
        for h in self.hists:
            if h=='reft':
                outfile.WriteObject(self.hists[h], self.hists[h].GetName(),'kOverwrite')

            elif len(h.split('_'))==3:
                for additionalkey in additionalkeys:
                    if h.find(additionalkey)==-1: continue
                    planekey, bar, key = h.split('_')
                    
                    hist=self.hists[h]
                    if not hasattr(outfile, additionalkey): folder=outfile.mkdir(additionalkey)
                    else: folder=outfile.Get(additionalkey)
                    folder.cd()
                    hist.Write(h, 2) # The 2 means it will overwrite a hist of the same name  

            elif len(h.split('_'))==2:
                key, plane = h.split('_')
                hist=self.hists[h]
                if not hasattr(outfile, 'xL-cdt'): folder=outfile.mkdir('xL-cdt')
                else: folder=outfile.Get('xL-cdt')
                folder.cd()
                hist.Write(h, 2) # The 2 means it will overwrite a hist of the same name                  

        outfile.Close()
        print(f'{len(self.hists)} histograms saved to {outfilename}')        

    # def yresidual3(self, detID, pos, mom):
    def yresidual3(self, detID):
        self.tw.MuFilter.GetPosition(detID,A,B)
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
            modaltime=hist.GetBinCenter(hist.GetMaximumBin())
            mean, stddev=hist.GetMean(), hist.GetStdDev()
            # print(f'td: {td}, mean +/- 2*std. dev.: {mean-2*stddev}, {mean+2*stddev}')
            if td < modaltime-2*stddev or td > modaltime+2*stddev: return False
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

    def reftcut(self,reft):
        if reft==-6237.5: return 0
        reftns=reft*self.TDC2ns
        if reftns<13.75 or reftns>15.45: return 0
        return reftns 
    
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
                      
