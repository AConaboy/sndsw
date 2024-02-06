#!/usr/bin/env python
import ROOT, csv, os
import numpy as np

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
        self.highQDCthreshold=60

        self.A, self.B = ROOT.TVector3(), ROOT.TVector3()
        
        if options.signalpartitions: self.Loadnumuevents()

        if tw.numuStudy: 
            self.numuStudy=True
            self.barycentres={}
        
    def Loadnumuevents(self):
        numusignalevent_filepath = '/afs/cern.ch/work/a/aconsnd/numusignalevents.csv'
        with open(numusignalevent_filepath, 'r') as f:
            reader=csv.reader(f)
            nu_mu_data=[r for r in reader]
        self.nu_mu_events={int(x[0]):(int(x[1]), int(x[2])) for x in nu_mu_data}        

    def FillHists(self, hits):
        
        """ 
        Hits passed here will be in events with a DS trackf
        I then check if they are in good agreement with the DS track in y.
        For this code I do NOT want hits associated with the track. 
        I will then check if they have atleast 4 (6) US (veto) SiPMs firing
        before agreeing to analyse the SiPM data
        
        For non-track correlated hits I will:
            1. Plot the fraction of fired SiPMs in the hit
            2. Find the next event where the each of the 4 small SiPMs fire. 
            3. Correlate the clock cycles with that delay with the number of fired channels on the PCB
        """

        self.nSiPMsPCB = self.muAna.GetFiredSiPMsOnPCBs(hits)
        
        for hit in hits: # Looking now at all hits, should start here incase d(phi)z between hadronic shower and the muon is sub 1 bar
            if not hit.isValid(): continue
            
            detID = hit.GetDetectorID()

            s, p, b= self.muAna.parseDetID(detID)
            if s==3:continue

            # Fill histogram with the fraction of the SiPMs that fire in this hit.
            fractionSiPMs_histname=f'frac_SiPMs'
            if not fractionSiPMs_histname in self.hists:
                title = 'Fraction of SiPMs that fire for hits in US;N_{fired} / N_{expected};Counts'
                self.hists[fractionSiPMs_histname] = ROOT.TH1F(fractionSiPMs_histname, title, 16, 0, 1)
            
            qdcs=hit.GetAllSignals()
            found_large_SiPMs=[]
            found_small_SiPMs=[]
            for i in qdcs:
                if i[0] in (2,5,10,13): found_small_SiPMs.append(i[0])
                else: found_large_SiPMs.append(i[0])
                
            if any([len(qdcs)==0, sum([i[1] for i in qdcs])==0 ]):
                print([(i[0],i[1]) for i in qdcs])
            
            frac_n = ( len(found_large_SiPMs) + len(found_small_SiPMs) ) / 16
            self.hists[fractionSiPMs_histname].Fill(frac_n) 

            if self.options.SmallSiPMcheck:
                """
                Soooo now I make the delta evt_timestamp histogram.
                Fill it with 0 for each small SiPM that fires in the hit.
                Call method to search next 100 events for the small SiPMs
                """
                
                next_smallSiPM_timestamp_histname=f'next_smallSiPM_timestamp'
                if not next_smallSiPM_timestamp_histname in self.hists:
                    title = 'Clock cycles until expected small SiPM next fires;#Delta clock cycles [n];Counts'
                    self.hists[next_smallSiPM_timestamp_histname] = ROOT.TH1I(next_smallSiPM_timestamp_histname, title, 50, 0, 50)
                
                # For the small SiPMs in this event, fill dt hist with 0 and remove it from the list of small SiPMs
                for smallSiPM in found_small_SiPMs: 
                    self.hists[next_smallSiPM_timestamp_histname].Fill(0)
                    
                """
                Look for the other small SiPMs over the next 100 events 
                Fill histograms of:
                    
                    1. Correlating the delay with nSiPMs/PCB 
                    2. Correlating the delay with QDC 
                """
                
                missing_small_SiPMs = [i for i in [2,5,10,13] if not i in found_small_SiPMs]
                
                if len(missing_small_SiPMs) != 0:
                    x=self.FindExpectedSmallSiPMs(detID, qdcs, missing_small_SiPMs)
                    if x==-999: continue # Skip event if no large SiPMs fire.
            
            ### Apply time-walk correction to all times in the hit.
            correctedtimes=self.muAna.GetCorrectedTimes(hit)
            times={'left':{}, 'right':{}}

            for i in correctedtimes: 
                SiPM, correctedtime = i
                fixed_ch=self.muAna.MakeFixedCh((s,p,b,SiPM))

                if not fixed_ch in self.muAna.alignmentparameters: continue
                d=self.muAna.alignmentparameters[fixed_ch]

                qdc=self.muAna.GetChannelVal(SiPM, qdcs)
                side=self.muAna.GetSide(f'{detID}_{SiPM}')
                times[side][SiPM]=self.tw.TDS0 - correctedtime - d[0]
                
                QDC_histname = f'US-SiPMQDC'
                if not QDC_histname in self.hists:
                    title='QDC_{SiPM} for large SiPMs in #nu_{#mu} candidates hits'
                    self.hists[QDC_histname]=ROOT.TH1F(QDC_histname,title, 110, -10, 100)
                self.hists[QDC_histname].Fill(qdc)

            if any([len(times[i])==0 for i in ('left', 'right')]): 
                print(f'No corrected and aligned times determined for this hit')
                continue
                # return
            averages={side:sum(times[side].values()) / len(times[side]) for side in times}
            averagetime = 1/2 * sum(averages.values())
            deltatime = averages['left'] - averages['right']
            
            ### Don't need to worry about nSiPMs because I have already cut on it
            baraveragecscint_left, baraveragecscint_right = self.muAna.GetBarAveragecscint(self.TWCorrectionRun, detID, 'corrected')

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

            # radiusestimate=f'DISradius_{self.subsystemdict[s]}'
            # if not radiusestimate in self.hists:
            #     title='Muon DIS shower radius estimate;L - c_{bar}(#frac{1}{t_{L}} - #frac{1}{t_{R}}) [cm];Counts'
            #     self.hists[radiusestimate]=ROOT.TH1F(radiusestimate, title, 200, -50, 50)
            # r=0.5*(self.barlengths[s] - baraveragecscint[0]*(1/averages['left'] - 1/averages['left']) )
            # self.hists[radiusestimate].Fill(r)    

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

    def FindExpectedSmallSiPMs(self, fDetID, qdcs, smallSiPMs):
        
        # Time of event in clock cycles since run began 
        
        # Use options.signalpartitions to get the right runNumber
        # Now this function will only work for nu_mu data
        
        large_SiPMQDCs=[i[1] for i in qdcs if not i[0] in (2,5,10,13)]
        if len(large_SiPMQDCs)==0: 
            print(f'Que??')
            return -999 # Weird that this happens?
        
        avg_large_SiPMQDC = sum(large_SiPMQDCs) / len(large_SiPMQDCs)

        n_missing = len(smallSiPMs)
        n_found = 0
        
        sndeventheader = self.tw.M.eventTree.EventHeader
        runNr = self.tw.M.eventTree.EventHeader.GetRunId()
        large_eventtime = sndeventheader.GetEventTime() # clock cycles since run started

        if self.options.signalpartitions: # If I'm looking at mu_nu candidates
            n_partition = list(self.options.signalpartitions.keys()).index(f'{str(runNr).zfill(6)}')
            event_in_partition = self.nu_mu_events[runNr][0] % 1e6 # self.nu_mu_events[runNr] stores: (event in partition, interaction wall)
            
            # signal_event = int(n_partition * 1e6 + self.nu_mu_events[runNr] % 1e6) 
            signal_event = int(n_partition * 1e6 + event_in_partition) 
        else: signal_event = self.tw.M.EventNumber
        
        # Loop over subsequent 100 events
        for evt in range(signal_event+1, signal_event+101):
            
            d_event = evt - signal_event

            self.tw.M.GetEvent(evt)
            hits=self.tw.M.eventTree.Digi_MuFilterHits
            
            # Check that the detID I need is in the event
            for h in hits: 
                detID = h.GetDetectorID() 
                
                # Check if required detID is in this event
                if detID != fDetID: continue 
                
                # Check if any required small SiPM is in this event
                SiPMs = [i[0] for i in h.GetAllSignals()]
                if not any([i in smallSiPMs for i in SiPMs]): continue 
                
                # Fill delta event time hist
                sndeventheader = self.tw.M.eventTree.EventHeader
                small_eventtime = sndeventheader.GetEventTime()
                d_clockcycles = small_eventtime - large_eventtime

                signals=h.GetAllSignals()
                for x in signals:
                    
                    SiPM, QDC = x 
                    if SiPM not in smallSiPMs: continue # If a required small SiPM doesn't fire in this hit then continue

                    smallSiPMs.remove(SiPM) # Don't double count small SiPM hits!
                    n_found += 1
                    
                    fixed_ch = f'{detID}_{SiPM}'
                    side = self.muAna.GetSide(fixed_ch)
                    
                    # Fill delta t hist for each small SiPM hit found in the delayed event
                    next_smallSiPM_timestamp_histname=f'next_smallSiPM_timestamp'
                    self.hists[next_smallSiPM_timestamp_histname].Fill(d_clockcycles)
                    
                    delayedSmallSiPM_dtvQDC = f'delayedSmallSiPM_{SiPM}_dtvQDC'
                    if not delayedSmallSiPM_dtvQDC in self.hists:
                        title='QDC v delta clock cycles for small SiPM '+str(SiPM)+' in '+str(detID)+';#Delta clock cycles (small SiPMs, #nu_{#mu} candidate);QDC [a.u]'
                        self.hists[delayedSmallSiPM_dtvQDC] = ROOT.TH2F(delayedSmallSiPM_dtvQDC, title, 50, 0, 50, 60, -10, 50)
                    self.hists[delayedSmallSiPM_dtvQDC].Fill(d_clockcycles, QDC) 
                    
                    s,p,b = self.muAna.parseDetID(detID)
                    SiPMsOnPCBkey = f'{p}_{side}'
                    SiPMsOnPCB = self.nSiPMsPCB[SiPMsOnPCBkey]
                    
                    delayedSmallSiPM_dtvnSiPMs = f'SiPMsOnPCB_dtvnSiPMs'
                    if not delayedSmallSiPM_dtvnSiPMs in self.hists:
                        title = 'Fired SiPMs on PCB v delay of small SiPM in clock cycles;#Delta clock cycles (small SiPMs, #nu_{#mu} candidate);Fired SiPMs on PCB'
                        self.hists[delayedSmallSiPM_dtvnSiPMs] = ROOT.TH2F(delayedSmallSiPM_dtvnSiPMs, title, 100, 0, 100, 50, 0, 50)
                    self.hists[delayedSmallSiPM_dtvnSiPMs].Fill(d_clockcycles, SiPMsOnPCB)
                    
                    delayedSmallSiPM_dtvnSiPMs = f'SiPMsOnPCB_dtvnSiPMs-highQDC'
                    if not delayedSmallSiPM_dtvnSiPMs in self.hists:
                        title = '#splitline{Fired SiPMs on PCB v delay of small SiPM in clock cycles}{average large SiPM QDC > '+str(self.highQDCthreshold)+'};#Delta clock cycles (small SiPMs, #nu_{#mu} candidate);Fired SiPMs on PCB'
                        self.hists[delayedSmallSiPM_dtvnSiPMs] = ROOT.TH2F(delayedSmallSiPM_dtvnSiPMs, title, 100, 0, 100, 50, 0, 50)
                    if avg_large_SiPMQDC>=self.highQDCthreshold: self.hists[delayedSmallSiPM_dtvnSiPMs].Fill(d_clockcycles, SiPMsOnPCB)                    
                    
                    largeSiPMQDC_dtvQDC = f'largeSiPMQDC_dtvQDC'
                    if not largeSiPMQDC_dtvQDC in self.hists:
                        title = 'Small SiPM delay in clock cycles v average QDC of the large SiPMs in the US hit;#Delta clock cycles (small SiPMs, #nu_{#mu} candidate);Average QDC per large SiPM [a.u]'
                        self.hists[largeSiPMQDC_dtvQDC] = ROOT.TH2F(largeSiPMQDC_dtvQDC, title, 100, 0, 100, 200, -20, 180) 
                    self.hists[largeSiPMQDC_dtvQDC].Fill(d_clockcycles, avg_large_SiPMQDC)
                    
                    largevsmall_QDCvQDC = f'largevsmall_QDCvQDC'
                    if not largevsmall_QDCvQDC in self.hists:
                        title = 'Correlation of average large SiPM QDC v small SiPM QDC for US hits in #nu_{#mu} candidates;Average QDC per large SiPM [a.u];Small SiPM QDC [a.u]'
                        self.hists[largevsmall_QDCvQDC] = ROOT.TH2F(largevsmall_QDCvQDC, title, 200, -20, 180, 60, -10, 50) 
                    self.hists[largevsmall_QDCvQDC].Fill(avg_large_SiPMQDC, QDC)
                    
                    smallQDCvdt = f'smallQDCvdt_dtvQDC'
                    if not smallQDCvdt in self.hists:
                        title = 'Correlation of small SiPM delay in clock cycles v small SiPM QDC for US hits in #nu_{#mu} candidates;#Delta clock cycles (small SiPMs, #nu_{#mu} candidate);Small SiPM QDC [a.u]'
                        self.hists[smallQDCvdt] = ROOT.TH2F(smallQDCvdt, title, 100, 0, 100, 60, -10, 50) 
                    self.hists[smallQDCvdt].Fill(d_clockcycles, QDC) 
                    
                    smallQDCvdt = f'smallQDCvdt_dtvQDC-highQDC'
                    if not smallQDCvdt in self.hists:
                        title = '#splitline{Correlation of small SiPM delay in clock cycles v small SiPM QDC for US hits in #nu_{#mu} candidates}{Average large SiPM QDC>'+str(self.highQDCthreshold)+'};#Delta clock cycles (small SiPMs, #nu_{#mu} candidate);Small SiPM QDC [a.u]'
                        self.hists[smallQDCvdt] = ROOT.TH2F(smallQDCvdt, title, 100, 0, 100, 60, -10, 50) 
                    if avg_large_SiPMQDC>=self.highQDCthreshold: self.hists[smallQDCvdt].Fill(d_clockcycles, QDC)                                        
          
        fractionMissingSmallFoundHighQDC = f'fractionMissingSmallFound-highQDC'
        if not fractionMissingSmallFoundHighQDC in self.hists:
            title = '#splitline{Fraction of small SiPMs missing from signal candidate hit}{which where found in the next 100 events. Average large SiPM QDC>'+str(self.highQDCthreshold)+'};Fraction exp. small SiPM firing within 100 clock cycles;Counts'
            self.hists[fractionMissingSmallFoundHighQDC] = ROOT.TH1F(fractionMissingSmallFoundHighQDC, title, 6, 0, 1.2)
        if avg_large_SiPMQDC>self.highQDCthreshold: self.hists[fractionMissingSmallFoundHighQDC].Fill(n_found/n_missing)

        fractionMissingSmallFound = f'fractionMissingSmallFound'
        if not fractionMissingSmallFound in self.hists:
            title = '#splitline{Fraction of small SiPMs missing from signal candidate hit}{which where found in the next 100 events};Fraction exp. small SiPM firing within 100 clock cycles;Counts'
            self.hists[fractionMissingSmallFound] = ROOT.TH1F(fractionMissingSmallFound, title, 6, 0, 1.2)
        self.hists[fractionMissingSmallFound].Fill(n_found/n_missing)        

        # Very important to load the signal event back!
        self.tw.M.GetEvent(signal_event)

    def ShowerDirection(self, hits):
        
        self.all_times={i:{} for i in range(5)}
        self.all_qdcs={i:{} for i in range(5)}
        
        # For numu study, 
        if self.numuStudy: 
            self.runId = self.tw.M.eventTree.EventHeader.GetRunId()
            if not int(self.runId) in self.barycentres: self.barycentres[self.runId] = {}

        """
        1. Plot shower baricentre in y (x) by using a QDC weighting, aligned timing
        2. First consider all hits 
        3. For each event, plot the baricentre in x and y at each plane
        """

        for hit in hits:
            if not hit.isValid(): continue

            detID = hit.GetDetectorID()
            s,p,b = self.muAna.parseDetID(detID)
            
            if s!=2: continue

            # Using statistical method: comparing the previously computed dy residual
            if self.tw.yresidual3(detID): self.trackrelated = True
            else: self.trackrelated=False

            # Extrapolate DS track to zPos for this plane, then simply not considering the bar which the DS track extrapolates to for the shower direction
            if self.muAna.GetExtrapolatedBarDetID(p) == detID: self.trackInBar = True
            else: self.trackInBar = False

            # run numusignals.py with --dycut to require US hits are NOT correlated in y with the DS track
            if self.options.dycut and self.trackrelated: continue
            elif self.options.notDSbar and self.trackInBar: continue

            self.DetermineAverages(hit)

        self.FillBarycentrePlots()

    def DetermineAverages(self, hit):
        
        detID = hit.GetDetectorID()
        s,p,b = self.muAna.parseDetID(detID)
        
        ### Apply time-walk correction to all times in the hit.
        qdcs = hit.GetAllSignals()

        # nSiPMs criteria on hits. Should replace with a histogram check
        if len(qdcs) < 11: return

        qdc_sides = {'left':[], 'right':[]}
        
        for idx in qdcs: 
            SiPM, qdc = idx 
            if SiPM < 8: qdc_sides['left'].append(qdc)
            else: qdc_sides['right'].append(qdc)

        correctedtimes=self.muAna.GetCorrectedTimes(hit, mode='aligned')
        
        aligned_times={'left':{}, 'right':{}}
        # for i in correctedtimes: 
        for i in aligned_times:
            SiPM, correctedtime = i
            if SiPM>8: aligned_times.append([SiPM, correctedtime])

        if any([len(aligned_times[i])==0 for i in ('left', 'right')]): return
        averages={side:sum(aligned_times[side].values()) / len(aligned_times[side]) for side in aligned_times}
        fastest = {side:min(aligned_times[side].values()) for side in aligned_times}
        averagetime = 1/2 * sum(averages.values())

        delta_averagetime = averages['left'] - averages['right']
        delta_fastesttime = fastest['left'] - fastest['right']
        
        self.all_times[p][detID] = averages
        self.all_qdcs[p][detID] = qdc_sides

    def FillBarycentrePlots(self):

        for plane in range(5):

            Nfiredbars = len(self.all_times[plane])
            times, qdcs = self.all_times[plane], self.all_qdcs[plane]

            # Extrapolate DS track to this plane
            zEx = self.zPos['MuFilter'][20+plane]
            lam = (zEx-self.tw.pos.z())/self.tw.mom.z()
            yEx = self.tw.pos.y() + lam*self.tw.mom.y()
            xEx = self.tw.pos.x() + lam*self.tw.mom.x()

            # Maintain dictionary stucture
            if any([len(data) == 0 for data in (times, qdcs)]): 
                if self.numuStudy: self.barycentres[self.runId][plane] = [[(-999,-999), (-999,-999)],[xEx, yEx]]
                continue

            # Both return -999 if fucked
            y_barycentre = self.GetYBarycentre(qdcs)
            x_barycentre = self.GetXBarycentre(qdcs, times, Nfiredbars)

            if self.numuStudy: self.barycentres[self.runId][plane] = [[x_barycentre, y_barycentre], [xEx, yEx]]

            deltay = yEx - y_barycentre[0]
            deltax = xEx - x_barycentre[0]

            # if self.numuStudy:
            #     self.barycentres[runId][plane] = [[x_barycentre, y_barycentre], []]
            
            y_barycentre_hist = f'y-barycentre_plane{plane}'
            if not y_barycentre_hist in self.hists:
                title = f'Barycentre in y from QDC weighting;y-barycentre [cm];Counts'
                self.hists[y_barycentre_hist] = ROOT.TH1F(y_barycentre_hist, title, 150, 0, 150)
            self.hists[y_barycentre_hist].Fill(y_barycentre[0])

            x_barycentre_hist = f'x-barycentre_plane{plane}'
            if not x_barycentre_hist in self.hists:
                title = 'Barycentre in x: 1/2*c_{scint}*(t_{L}+t_{R});x-barycentre [cm];Counts'
                self.hists[x_barycentre_hist] = ROOT.TH1F(x_barycentre_hist, title, 200, -100, 100)
            self.hists[x_barycentre_hist].Fill(x_barycentre[0])            

            deltaymuonbarycentre = f'deltay-muonbarycentre_plane{plane}'
            if not deltaymuonbarycentre in self.hists:
                title='Difference in y between muon and barycentre from QDC for plane '+str(plane+1)+';#Delta y (DS track, QDC barycentre);Counts'
                self.hists[deltaymuonbarycentre] = ROOT.TH1F(deltaymuonbarycentre, title, 140, -70, 70)
            self.hists[deltaymuonbarycentre].Fill(deltay)
        
            deltaxmuonbarycentre = f'deltax-muonbarycentre_plane{plane}'
            if not deltaxmuonbarycentre in self.hists:
                title='Difference in x between muon and barycentre from timing for plane '+str(plane+1)+';#Delta x (DS track, QDC barycentre);Counts'
                self.hists[deltaxmuonbarycentre] = ROOT.TH1F(deltaxmuonbarycentre, title, 140, -70, 70)
            self.hists[deltaxmuonbarycentre].Fill(deltax)

    def GetYBarycentre(self, qdcs):
        
        barQDCs = {k:0 for k in qdcs}
        for detID in qdcs:
            barQDCs[detID] = sum( [sum(innerlist) for innerlist in qdcs[detID].values()] )

        sumQDC = sum([i for i in barQDCs.values()])

        A,B = ROOT.TVector3(), ROOT.TVector3()
        weighted_y_positions = []
        for key in qdcs:
            self.MuFilter.GetPosition(key, A, B)
            weighted_y_positions.append( barQDCs[key]/sumQDC * 1/2 * (A.y() +B.y()) )

        if not len(weighted_y_positions)==0: return sum(weighted_y_positions)/len(weighted_y_positions), 0
        else: return -999, -999
        
    def GetXBarycentre(self, qdcs, times, Nfiredbars):
        
        """
        Determine x-barycentre using bar with largest QDC
        Already, data contains the qdcs of the larger side max( sum(right_qdcs), sum(left_qdcs) )
        So I can use detID:[qdcs] to find the bar with the largest QDC deposit
        """
        x = {}
        for detID in times:
            s,p,b = self.muAna.parseDetID(detID)

            averages = times[detID]       
            if len(averages) == 0: return -999,-999

            average_left_cscint, average_right_cscint = self.muAna.GetBarAveragecscint(self.runNr, detID, 'corrected')

            if detID == self.muAna.GetExtrapolatedBarDetID(p): 
                xEx, yEy, zEx = self.muAna.GetExtrapolatedPosition(p)
                self.MuFilter.GetPosition(detID, self.A, self.B) # A is left, B is right
                shower_side = 'right' if abs(self.A.x() - xEx)>abs(self.B.x() - xEx) else 'left'
                c_showerside = average_left_cscint if shower_side=='left' else average_right_cscint
                x_shower = averages[shower_side] * c_showerside[0]
                x_barycentre = x_shower + Nfiredbars/2*6 # 6 cm is from MuFilter.GetConfParF('MuFilter/UpstreamBarY')
                x_barycentre_uncertainty = 3
                x[detID] = [x_barycentre, x_barycentre_uncertainty]
            else:
            
                x_barycentre = 0.5*(average_left_cscint[0]+average_right_cscint[0]) * 0.5 * (averages['right'] + averages['left'])

                # Approximate uncertainty right now as:  [ (dcscint/cscint)**2 + 1/2*(150ps/(tL+tR)**2 )]**(1/2)
                # rel_cscint_uncetainty_sq = 1/4 *( (average_left_cscint[1]/average_left_cscint[0])**2 + (average_right_cscint[1]/average_right_cscint[0])**2 ) / (average_left_cscint[0] + average_right_cscint[0])
                rel_cscint_uncetainty_sq = 1/4 *( (average_left_cscint[1]/average_left_cscint[0])**2 + (average_right_cscint[1]/average_right_cscint[0])**2 )
                
                rel_tL_tR_uncetainty_sq = 1/np.sqrt(2) * 0.150 **2

                x_barycentre_uncertainty_sq = x_barycentre**2 * (rel_cscint_uncetainty_sq + rel_tL_tR_uncetainty_sq)
                x_barycentre_uncertainty = np.sqrt(x_barycentre_uncertainty_sq)

                x[detID] = [x_barycentre, x_barycentre_uncertainty]

        mean_value = sum(item[0] for item in x.values()) / len(x)
        mean_uncertainty = sum(item[1] for item in x.values()) / len(x)
        return mean_value, mean_uncertainty
        
    def GetAlignedTimes(self, detID, correctedtimes, qdcs):
        times={'left':{}, 'right':{}}
        for i in correctedtimes: 
            SiPM, correctedtime = i
            # fixed_ch=self.muAna.MakeFixedCh((s,p,b,SiPM))
            fixed_ch=f'{detID}_{SiPM}'

            if not fixed_ch in self.muAna.alignmentparameters: continue
            d=self.muAna.alignmentparameters[fixed_ch]

            qdc=self.muAna.GetChannelVal(SiPM, qdcs)
            side=self.muAna.GetSide(f'{detID}_{SiPM}')
            times[side][SiPM]=self.tw.TDS0 - (correctedtime - d[0])
        return times
        
    def WriteOutHistograms(self):

        outfilename=f'{self.outpath}splitfiles/run{self.runNr}/ShowerProfiles/ShowerProfiles_{self.options.nStart}.root'
        if os.path.exists(outfilename): outfile=ROOT.TFile.Open(outfilename, 'recreate')
        else: outfile=ROOT.TFile.Open(outfilename, 'create')

        additionalkeys=['averagetime', 'deltatime', 'sidetime', 'AlignedSiPMtime', 'timingxt', 'y-barycentre', 'x-barycentre', 'deltay-muonbarycentre', 'deltax-muonbarycentre']
        for h in self.hists:
            if h=='TDS0':
                outfile.WriteObject(self.hists[h], self.hists[h].GetName(),'kOverwrite')

            elif len(h.split('_'))>1:
                for additionalkey in additionalkeys:
                    if h.find(additionalkey)==-1: continue
                    
                    hist=self.hists[h]
                    if not hasattr(outfile, additionalkey): folder=outfile.mkdir(additionalkey)
                    else: folder=outfile.Get(additionalkey)
                    folder.cd()
                    hist.Write(h, 2) # The 2 means it will overwrite a hist of the same name

        outfile.Close()
        print(f'{len(self.hists)} histograms saved to {outfilename}')    

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
                      
