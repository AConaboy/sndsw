#!/usr/bin/env python
import ROOT

class ShowerProfiles(object):

    def __init__(self, options, tw):
       
        # self.state='corrected'
        # options.state='corrected'
        self.options=options
        self.tw=tw
        self.runNr = tw.runNr 
        self.muAna = tw.muAna
        self.timealignment=tw.timealignment
        self.TWCorrectionRun=tw.TWCorrectionRun

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

        self.MuFilter = tw.MuFilter
        self.barlengths = self.muAna.BuildBarLengths(self.MuFilter)

        self.freq=160.316E6
        self.TDC2ns=1E9/self.freq

        self.xrefs=tw.xrefs
        self.sides=('left', 'right')

        self.hists=tw.hists
        
        self.sigmatds0=0.263 # ns 

    def FillHists(self, hits):
        """ 
        Hits passed here will be in events with a DS track
        I then check if they are in good agreement with the DS track in y.
        For this code I do NOT want hits associated with the track. 
        I will then check if they have atleast 4 (6) US (veto) SiPMs firing
        before agreeing to analyse the SiPM data
        """
        nontrackhits=[]
        for hit in hits:

            detID = hit.GetDetectorID()
            s,p,b = self.muAna.parseDetID(detID)
            if s==3: continue
            
            if self.tw.yresidual3(detID): continue
            
            # Currently the SiPM cut doesn't account for the few detIDs with missing channels, 
            # so I probably have a bar where I won't see data for 1 bar in the veto and maybe 1 in the US
            nLeft, nRight=self.muAna.GetnFiredSiPMs(hit)
            if not self.nSiPMscut(hit, nLeft, nRight): continue 
            
            nontrackhits.append(hit)

        nt_multiplicity='ExtraHitsMultiplicity'
        if not nt_multiplicity in self.hists:
            nt_multiplicitytitle = f'Multiplicity of non-track related hits in DS track events;Multiplicity;Counts'
            self.hists[nt_multiplicity]=ROOT.TH1F(nt_multiplicity, nt_multiplicitytitle, 6, 0, 6)
        self.hists[nt_multiplicity].Fill(len(nontrackhits))

        for nt_hit in nontrackhits:

            detID=nt_hit.GetDetectorID()
            s, p, b= self.muAna.parseDetID(detID)

            ### Apply time-walk correction to all times in the hit.
            ### GetCorrectedTimes(muAna, hit, x=0, mode='unaligned')

            correctedtimes=self.muAna.GetCorrectedTimes(nt_hit)
            qdcs = nt_hit.GetAllSignals()
            # print(f'detID: {detID}, len corrected times: {len(correctedtimes)}')
            times={'left':{}, 'right':{}}

            for i in correctedtimes: 
                SiPM, correctedtime = i
                fixed_ch=self.muAna.MakeFixedCh((s,p,b,SiPM))

                if not fixed_ch in self.muAna.alignmentparameters: continue
                d=self.muAna.alignmentparameters[fixed_ch]

                qdc=self.muAna.GetChannelVal(SiPM, qdcs)
                side=self.muAna.GetSide(f'{detID}_{SiPM}')
                times[side][SiPM]=self.tw.TDS0 - correctedtime - d[0]

            if any([len(times[i])==0 for i in ('left', 'right')]): return
            averages={side:sum(times[side].values()) / len(times[side]) for side in times}
            averagetime = 1/2 * sum(averages.values())
            deltatime = averages['left'] - averages['right']
            
            ### Don't need to worry about nSiPMs because I have already cut on it
            baraveragecscint = self.muAna.GetBarAveragecscint(self.TWCorrectionRun, detID, 'corrected')

            averagebartimehistname=f'averagetime_{detID}_aligned'
            if not averagebartimehistname in self.hists:
                title=self.subsystemdict[s]+' plane '+str(p+1)+' bar '+str(b+1)+' average of aligned times;#frac{1}{2}#times(t^{tw corr}_{left}+t^{tw corr}_{right}) [ns];Counts'
                self.hists[averagebartimehistname]=ROOT.TH1F(averagebartimehistname, title, 2000, -5, 5)
            self.hists[averagebartimehistname].Fill(averagetime)
            
            deltabartimehistname=f'deltatime_{detID}_aligned'
            if not deltabartimehistname in self.hists:
                title=self.subsystemdict[s]+' plane '+str(p+1)+' bar '+str(b+1)+' difference of average aligned time of each side;t^{tw corr}_{left} - t^{tw corr}_{right} [ns];Counts'
                self.hists[deltabartimehistname]=ROOT.TH1F(deltabartimehistname, title, 2000, -5, 5)
            self.hists[deltabartimehistname].Fill(deltatime)

            for side in ('left', 'right'):
                averagebarsidetimehistname=f'sidetime_{detID}-{side}_aligned'
                if not averagebarsidetimehistname in self.hists:
                    title=self.subsystemdict[s]+' plane '+str(p+1)+' bar '+str(b+1)+' average aligned time from '+side+' side;'+ side+' side average of t^{tw corr}_{'+side+'} [ns];Counts'
                    self.hists[averagebarsidetimehistname]=ROOT.TH1F(averagebarsidetimehistname, title, 2000, -5, 5)
                self.hists[averagebarsidetimehistname].Fill(averages[side])

            radiusestimate=f'DISradius_{self.subsystemdict[s]}'
            if not radiusestimate in self.hists:
                title='Muon DIS shower radius estimate;L - c_{bar}(#frac{1}{t_{L}} - #frac{1}{t_{R}}) [cm];Counts'
                self.hists[radiusestimate]=ROOT.TH1F(radiusestimate, title, 200, -50, 50)
            r=0.5*(self.barlengths[s] - baraveragecscint[0]*(1/averages['left'] - 1/averages['left']) )
            self.hists[radiusestimate].Fill(r)

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
                      
