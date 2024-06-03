#!/usr/bin/env python
import ROOT, os

class SystemAlignment(object):

    def __init__(self, options, tw):
       
        self.options=options
        self.tw=tw
        self.simulation = tw.simulation
        # if self.tw.mode=='reconstructmuonposition': 
        self.USbarlength = self.tw.MuFilter.GetConfParF('MuFilter/UpstreamBarX')
        self.runNr = tw.runNr 
        self.muAna = tw.muAna
        self.timealignment=tw.timealignment
        self.A, self.B = ROOT.TVector3(), ROOT.TVector3()

        self.afswork=tw.afswork
        if not self.simulation: self.outpath=tw.outpath
        else: self.outpath = options.path

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

        correctedtimes=self.muAna.GetCorrectedTimes(hit, self.tw.Ex.x(), mode='aligned')

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
                self.hists[SiPMtime]=ROOT.TH1F(SiPMtime,fulltitle, 400, -10, 10)

            self.hists[SiPMtime].Fill(time)

    def FillBarHists(self, hit):
        
        detID = hit.GetDetectorID()
        s, p, b = self.muAna.parseDetID(detID)
        self.tw.MuFilter.GetPosition(detID, self.A, self.B)
        x_midpoint = 0.5 * (self.A.x() + self.B.x())        

        correctedtimes = self.muAna.GetCorrectedTimes(hit, x=0, mode='aligned')
        tofcorrectedtimes = self.muAna.GetCorrectedTimes(hit, x=self.tw.Ex.x(), mode='aligned')
        print(f'Len atimes, tof atimes: {len(correctedtimes)}, {len(tofcorrectedtimes)}')
        aligned_times={'left':{}, 'right':{}}
        tofaligned_times={'left':{}, 'right':{}}

        for SiPM, alignedtime in correctedtimes:
            fixed_ch=self.muAna.MakeFixedCh((s,p,b,SiPM))
            side=self.muAna.GetSide(f'{detID}_{SiPM}')
            aligned_times[side][SiPM] = alignedtime

        for SiPM, tofatime in tofcorrectedtimes: 
            fixed_ch=self.muAna.MakeFixedCh((s,p,b,SiPM))
            side=self.muAna.GetSide(f'{detID}_{SiPM}')            
            tofaligned_times[side][SiPM] = tofatime

        if any([len(aligned_times[i])==0 for i in ('left', 'right')]): return
        if any([len(tofaligned_times[i])==0 for i in ('left', 'right')]): return
        averages = {side:(sum(aligned_times[side].values()) / len(aligned_times[side])) for side in aligned_times}
        fastest = {side:min(aligned_times[side].values()) for side in aligned_times}
        averagetime = 1/2 * sum(averages.values())

        cscint_left, cscint_right = self.muAna.GetBarAveragecscint(self.runNr, detID, 'corrected')
        delta_averagetime = averages['right'] - averages['left']
        delta_fastesttime = fastest['right'] - fastest['left']
        
        # For determining bar side-wise time resolution
        tofa_averages = {side:(sum(tofaligned_times[side].values()) / len(tofaligned_times[side])) for side in tofaligned_times}

        x='yAgreement' if self.tw.trackrelated else 'notyAgreement'

        # Looking at differences depending on whether the average or fastest times are used
        for timer in ('fastest', 'average'):
            dt = f'dt-{x}_plane{p}_{timer}'
            if not dt in self.hists:
                title='#Delta t{(left, right);#Delta t^{TW,DS}_{L}, t^{TW,DS}_{R} [ns];Counts'
                self.hists[dt]=ROOT.TH1F(dt,title,50, -5, 5)
            
            barycentrex_hist = f'xbarycentre-{x}_plane{p}_{timer}'
            if not barycentrex_hist in self.hists:
                title='Barycentre in x, #frac{1}{2}(x_{right}+x_{left});#frac{1}{2}(x_{right}+x_{left}) [cm];Counts'
                self.hists[barycentrex_hist] = ROOT.TH1F(barycentrex_hist, title, 110, -100, 10)

            lamx_hist = f'lambda-{x}_plane{p}_{timer}'
            if not lamx_hist in self.hists:
                title='Deposition width #{1}{2}(x_{left} + x_{right});#{1}{2}(x_{left} + x_{right}) [cm];Counts'
                self.hists[lamx_hist] = ROOT.TH1F(lamx_hist, title, 160, -90, 60)

            barycentrex_hist = f'xbarycentre-{x}_plane{p}_{timer}'
            if not barycentrex_hist in self.hists:
                title='Barycentre in x, #frac{1}{2}(x_{right}+x_{left});#frac{1}{2}(x_{right}+x_{left}) [cm];Counts'
                self.hists[barycentrex_hist] = ROOT.TH1F(barycentrex_hist, title, 110, -100, 10)

            if timer=='fastest':        
                self.hists[dt].Fill(delta_fastesttime)

                xL, xR = (x_midpoint-fastest['left']*cscint_left[0]), (fastest['right']*cscint_right[0]+x_midpoint) # Adding the x_midpoint translates into the physics FoR
                x_barycentre = 0.5*(xL + xR)
                lamx = 0.5*(xL - xR)

                self.hists[lamx_hist].Fill(lamx)
                self.hists[barycentrex_hist].Fill(x_barycentre)
            
            elif timer=='average':
                self.hists[dt].Fill(delta_averagetime)
                
                xL, xR = (x_midpoint-averages['left']*cscint_left[0]), (averages['right']*cscint_right[0]+x_midpoint) # Adding the x_midpoint translates into the physics FoR
                x_barycentre = 0.5*(xL + xR)
                lamx = 0.5*(xL - xR)

                self.hists[lamx_hist].Fill(lamx)
                self.hists[barycentrex_hist].Fill(x_barycentre)
        
        averagebartimehistname=f'averagetime_{detID}_aligned'
        if not averagebartimehistname in self.hists:
            title=self.subsystemdict[s]+' plane '+str(p+1)+' bar '+str(b+1)+' average of aligned times from each side as a function of x_{predicted};x_{predicted} [cm];'+'#frac{1}{2}#times(t^{DS,TW}_{left}+t^{DS,TW}_{right}) [ns];Counts'
            self.hists[averagebartimehistname]=ROOT.TH2F(averagebartimehistname, title, 110, -100, 10, 2000, -5, 5)
        self.hists[averagebartimehistname].Fill(self.tw.Ex.x(), averagetime)
        
        deltabartimehistname=f'deltatime_{detID}_aligned'
        if not deltabartimehistname in self.hists:
            title=self.subsystemdict[s]+' plane '+str(p+1)+' bar '+str(b+1)+' difference of average aligned time from each side as a function of x_{predicted};x_{predicted} [cm];'+'t^{DS,TW}_{left} - t^{DS,TW}_{right} [ns];Counts'
            self.hists[deltabartimehistname]=ROOT.TH2F(deltabartimehistname, title, 110, -100, 10, 2000, -5, 5)
        self.hists[deltabartimehistname].Fill(self.tw.Ex.x(), delta_fastesttime)

        for side in ('left', 'right'):
            averagebarsidetimehistname=f'sidetime_{detID}-{side}_aligned'
            if not averagebarsidetimehistname in self.hists:
                title=self.subsystemdict[s]+' plane '+str(p+1)+' bar '+str(b+1)+' average aligned + ToF corr time from '+side+' side as a function of x_{predicted};x_{predicted} [cm];'+ side+' side average of t^{DS,TW}_{'+side+'} [ns];Counts'
                self.hists[averagebarsidetimehistname]=ROOT.TH2F(averagebarsidetimehistname, title, 110, -100, 10, 2000, -5, 5)
            self.hists[averagebarsidetimehistname].Fill(self.tw.Ex.x(), tofa_averages[side])

    def ReconstructMuonPosition(self, hits):
        # self.MakeReconstructionHists()

        # planewise_data = {i:{} for i in range(5)}
        planewise_data = {}
        
        # Group hits by plane
        for hit in hits:
            detID = hit.GetDetectorID()
            s,p,b = self.muAna.parseDetID(detID)
            if not s==2 or not hit.isValid(): continue

            # I set tw.reft to the scifi average time before calling this method
            alignedtimes=self.muAna.GetCorrectedTimes(hit, mode='aligned')

            atimes_left = [i[1] for i in alignedtimes if i[0]<8]
            atimes_right = [i[1] for i in alignedtimes if i[0]>=8] 

            # Check nSiPMs left and right. Skip hit if less than 4 SiPMs on either left or right fire
            if len(atimes_left)<4 or len(atimes_right)<4:continue

            atimes_left_mean = sum(atimes_left)/len(atimes_left)
            atimes_right_mean = sum(atimes_right)/len(atimes_right)
            
            # Skip hit if abs( atimes_left(right) ) > L/2 / cscint_L(R)
            averagecscint_left, averagecscint_right = self.muAna.GetBarAveragecscint(self.tw.TWCorrectionRun, detID, 'corrected')
            if abs(atimes_left_mean) > self.USbarlength/2 / averagecscint_left[0] or abs(atimes_right_mean) > self.USbarlength/2 / averagecscint_right[0]: continue
            
            if not p in planewise_data: planewise_data[p]={}
            planewise_data[p][detID] = {}
            
            planewise_data[p][detID]['bar-QDC'] = self.muAna.GetTotalQDC(hit.GetAllSignals())

            planewise_data[p][detID]['atimes-left'] = atimes_left_mean
            planewise_data[p][detID]['atimes-right'] = atimes_right_mean

            planewise_data[p][detID]['cscint-left'] = averagecscint_left[0]
            planewise_data[p][detID]['cscint-right'] = averagecscint_right[0]

        for plane in planewise_data:
            pdata = planewise_data[plane]

            # Determine quantities for the bars now
            planeQDC = sum([ pdata[detID]['bar-QDC'] for detID in pdata ])

            weighted_ys = []
            x_barycentres = {}

            xEx, yEx, zEx = self.muAna.GetExtrapolatedPosition(plane)

            for detID in pdata:
                # Get weighted y-position for each hit that passes selection in this plane
                s,p,b = self.muAna.parseDetID(detID)
                
                if self.muAna.GetExtrapolatedBarDetID(p) == detID: trackInBar = True
                else: trackInBar=False
                
                barQDC = pdata[detID]['bar-QDC']
                self.tw.MuFilter.GetPosition(detID, self.A, self.B)
                y_pos = barQDC/planeQDC * 0.5 * (self.A.y() + self.B.y())
                x_midpoint = 0.5 * (self.A.x() + self.B.x())
                weighted_ys.append(y_pos) # Get weighted y-pos

                # Get x-barycentre
                cscint_left, cscint_right = pdata[detID]['cscint-left'],pdata[detID]['cscint-right']
                atimes_left, atimes_right = pdata[detID]['atimes-left'],pdata[detID]['atimes-right']

                # xL, xR = (x_midpoint-atimes_left*cscint_left), (atimes_right*cscint_right+x_midpoint) # Adding the x_midpoint translates into the physics FoR
                xL, xR = (x_midpoint+atimes_left*cscint_left), (-atimes_right*cscint_right+x_midpoint) # Adding the x_midpoint translates into the physics FoR
                x_barycentre = 0.5*(xL + xR)
                x_barycentres[detID] = x_barycentre

                for idx,xm in enumerate( (xL, xR, x_barycentre) ):

                    if idx==0:
                        deltabar_histname = f'{detID}-deltaxL'
                        deltabar_title = '#splitline{x-left - track x-position}{'+self.muAna.MakeHumanReadableDetID(detID)+'};#Delta (x_{left}, x_{Ex}) [cm];Counts'
                        bar_histname = f'{detID}-xL'
                        bar_title = '#splitline{x-left}{'+self.muAna.MakeHumanReadableDetID(detID)+'};x_{left} [cm];Counts'

                        deltaplane_histname = f'plane{p}-deltaxL'
                        deltaplane_title = '#splitline{x-left - track x-position}{plane '+str(p+1)+'};#Delta (x_{left}, x_{Ex}) [cm];Counts'

                        deltacorrplane_histname = f'plane{p}-dxLvxEx'
                        deltacorrplane_title = '#splitline{(x-left - track x-position) v track x-position}{plane '+str(p+1)+'};x_{track} [cm];#Delta (x_{left}, x_{Ex}) [cm]'

                        plane_histname = f'plane{p}-xL'
                        plane_title = '#splitline{x-left}{plane '+str(p+1)+'};x_{left} [cm];Counts' 
                        
                        plane_corrhistname = f'plane{p}-xLvxEx'
                        plane_corrtitle = '#splitline{x-left v track-x}{plane '+str(p+1)+'};x_{left} [cm];x_{track} [cm]'

                        deltatotal_histname = f'total-deltaxL'
                        deltatotal_title = 'x-left - track x-position for all bars in HCAL;#Delta (x_{left}, x_{Ex}) [cm];Counts'

                        deltacorrtotal_histname = f'total-dxLvxEx'
                        deltacorrtotal_title = '(x-left - track x-position) v track x-position for all bars in HCAL;x_{track} [cm];#Delta (x_{left}, x_{Ex}) [cm]'

                        total_histname = f'total-xL'
                        total_title = 'x-left for all bars in HCAL;x_{left} [cm];Counts'

                        total_corrhistname = f'total-xLvxEx'
                        total_corrtitle = 'x-left v track-x for all bars in HCAL;x_{left} [cm];x_{track} [cm]'

                    elif idx==1:
                        deltabar_histname = f'{detID}-deltaxR'
                        deltabar_title = '#splitline{x-right - track x-position}{'+self.muAna.MakeHumanReadableDetID(detID)+'};#Delta (x_{right}, x_{Ex}) [cm];Counts'
                        
                        bar_histname = f'{detID}-xR'
                        bar_title = '#splitline{x-right}{'+self.muAna.MakeHumanReadableDetID(detID)+'};x_{right} [cm];Counts'
                        
                        deltaplane_histname = f'plane{p}-deltaxR'
                        deltaplane_title = '#splitline{x-right - track x-position}{plane '+str(p+1)+'};#Delta (x_{right}, x_{Ex}) [cm];Counts'

                        deltacorrplane_histname = f'plane{p}-dxRvxEx'
                        deltacorrplane_title = '#splitline{(x-right - track x-position) v track x-position}{plane '+str(p+1)+'};x_{track} [cm];#Delta (x_{right}, x_{Ex}) [cm]'
                        
                        plane_histname = f'plane{p}-xR'
                        plane_title = '#splitline{x-right}{plane '+str(p+1)+'};x_{right} [cm];Counts'

                        plane_corrhistname = f'plane{p}-xRvxEx'
                        plane_corrtitle = '#splitline{x-right v track-x}{plane '+str(p+1)+'};x_{right} [cm];x_{track} [cm]'
                        
                        deltatotal_histname = f'total-deltaxR'
                        deltatotal_title = 'x-right - track x-position for all bars in HCAL;#Delta (x_{right}, x_{Ex}) [cm];Counts'

                        deltacorrtotal_histname = f'total-dxRvxEx'
                        deltacorrtotal_title = '(x-right - track x-position) v track x-position for all bars in HCAL;x_{track} [cm];#Delta (x_{right}, x_{Ex}) [cm]'

                        total_histname = f'total-xR'
                        total_title = 'x-right for all bars in HCAL;x_{right} [cm];Counts'

                        total_corrhistname = f'total-xRvxEx'
                        total_corrtitle = 'x-right v track-x for all bars in HCAL;x_{right} [cm];x_{track} [cm]'
                    
                    elif idx==2:
                        deltabar_histname = f'{detID}-deltaxB'
                        deltabar_title = '#splitline{x-barycentre - track x-position}{'+self.muAna.MakeHumanReadableDetID(detID)+'};#Delta (x_{barycentre}, x_{Ex}) [cm];Counts'
                        
                        bar_histname = f'{detID}-xB'
                        bar_title = '#splitline{x-barycentre}{'+self.muAna.MakeHumanReadableDetID(detID)+'};x_{barycentre} [cm];Counts'
                        
                        deltaplane_histname = f'plane{p}-deltaxB'
                        deltaplane_title = '#splitline{x-barycentre - track x-position}{plane '+str(p+1)+'};#Delta (x_{barycentre}, x_{Ex}) [cm];Counts'
                        
                        deltacorrplane_histname = f'plane{p}-dxBvxEx'
                        deltacorrplane_title = '#splitline{#Delta(x-barycentre - track x-position) v track x-position}{plane '+str(p+1)+'};x_{track} [cm];#Delta (x_{barycentre}, x_{Ex}) [cm]'

                        plane_histname = f'plane{p}-xB'
                        plane_title = '#splitline{x-barycentre}{plane '+str(p+1)+'};x_{barycentre} [cm];Counts'      

                        plane_corrhistname = f'plane{p}-xBvxEx'
                        plane_corrtitle = '#splitline{x-barycentre v track-x}{plane '+str(p+1)+'};x_{barycentre} [cm];x_{track} [cm]'
                        
                        deltatotal_histname = f'total-deltaxB'
                        deltatotal_title = 'x-barycentre - track x-position for all bars in HCAL;#Delta (x_{barycentre}, x_{Ex}) [cm];Counts'

                        deltacorrtotal_histname = f'total-dxBvxEx'
                        deltacorrtotal_title = '#Delta(x-barycentre - track x-position) v track x-position for all bars in HCAL;x_{track} [cm];#Delta (x_{barycentre}, x_{Ex}) [cm]'

                        total_histname = f'total-xB'
                        total_title = 'x-barycentre for all bars in HCAL;x_{barycentre} [cm];Counts'

                        total_corrhistname = f'total-xBvxEx'
                        total_corrtitle = 'x-barycentre v track-x for all bars in HCAL;x_{barycentre} [cm];x_{track} [cm]'
                    
                    # Fill dx hist per bar
                    if not bar_histname in self.hists:
                        self.hists[bar_histname] = ROOT.TH1F(bar_histname, bar_title, 120, -110, 10)
                    self.hists[bar_histname].Fill(xm)
                    # Fill dx hist per plane
                    if not plane_histname in self.hists:
                        self.hists[plane_histname] = ROOT.TH1F(plane_histname, plane_title, 120, -110, 10)
                    self.hists[plane_histname].Fill(xm) 
                    # Fill dx hist in system
                    if not total_histname in self.hists:
                        self.hists[total_histname] = ROOT.TH1F(total_histname, total_title, 120, -110, 10)
                    self.hists[total_histname].Fill(xm)

                    # Fill corr plots                 
                    if not plane_corrhistname in self.hists:
                        self.hists[plane_corrhistname] = ROOT.TH2F(plane_corrhistname, plane_corrtitle, 120, -110, 10, 120, -110, 10)
                    self.hists[plane_corrhistname].Fill(xm, xEx)
                    if not total_corrhistname in self.hists:
                        self.hists[total_corrhistname] = ROOT.TH2F(total_corrhistname, total_corrtitle, 120, -110, 10, 120, -110, 10)
                    self.hists[total_corrhistname].Fill(xm, xEx)

                    # Fill corr plots res v x                 
                    if not deltacorrplane_histname in self.hists:
                        self.hists[deltacorrplane_histname] = ROOT.TH2F(deltacorrplane_histname, deltacorrplane_title, 120, -110, 10, 100, -50, 50)
                    self.hists[deltacorrplane_histname].Fill(xEx,xm-xEx)
                    if not deltacorrtotal_histname in self.hists:
                        self.hists[deltacorrtotal_histname] = ROOT.TH2F(deltacorrtotal_histname, deltacorrtotal_title, 120, -110, 10, 100, -50, 50)
                    self.hists[deltacorrtotal_histname].Fill(xEx,xm-xEx)
                    
                    # Fill dx hist per bar
                    if not deltabar_histname in self.hists:
                        self.hists[deltabar_histname] = ROOT.TH1F(deltabar_histname, deltabar_title, 100, -50, 50)
                    self.hists[deltabar_histname].Fill(xm - xEx)
                    # Fill dx hist per plane
                    if not deltaplane_histname in self.hists:
                        self.hists[deltaplane_histname] = ROOT.TH1F(deltaplane_histname, deltaplane_title, 100, -50, 50)
                    self.hists[deltaplane_histname].Fill(xm - xEx) 
                    # Fill dx hist in system
                    if not deltatotal_histname in self.hists:
                        self.hists[deltatotal_histname] = ROOT.TH1F(deltatotal_histname, deltatotal_title, 100, -50, 50)
                    self.hists[deltatotal_histname].Fill(xm - xEx)  

                # Fill xL v xR in system
                xLvxRtotal_histname = f'total-xLvxR'
                if not xLvxRtotal_histname in self.hists:
                    title = f'total x-left v x-right;x-left [cm];x-right [cm];Counts'
                    self.hists[xLvxRtotal_histname] = ROOT.TH2F(xLvxRtotal_histname, title, 110, -100, 10, 110, -100, 10)
                self.hists[xLvxRtotal_histname].Fill(xL, xR)
                # Fill xL v xR in plane
                xLvxRplane_histname = f'plane{plane}-xLvxR'
                if not xLvxRplane_histname in self.hists:
                    title = '#splitline{x-left v x-right}{plane '+str(plane+1)+'};x-left [cm];x-right [cm];Counts'
                    self.hists[xLvxRplane_histname] = ROOT.TH2F(xLvxRplane_histname, title, 110, -100, 10, 110, -100, 10)
                self.hists[xLvxRplane_histname].Fill(xL, xR)                                  

            y_barycentre = sum(weighted_ys)      # Get plane y-barycentre
            if plane==0: print(pdata)

            # Fill yb hist in system
            total_histname = f'total-y'
            if not total_histname in self.hists:
                title = 'total US y-barycentre;y_{B} [cm];Counts'
                self.hists[total_histname] = ROOT.TH1F(total_histname, title, 100, -10, 90)
            self.hists[total_histname].Fill(y_barycentre)
            # Fill yb hist per plane
            plane_histname = f'plane{plane}-y'
            if not plane_histname in self.hists:
                title = '#splitline{y-barycentre}{plane '+str(plane+1)+'};y_{B} [cm];Counts'
                self.hists[plane_histname] = ROOT.TH1F(plane_histname, title, 100, -10, 90)
            self.hists[plane_histname].Fill(y_barycentre)

            # Fill dy hist in system
            deltatotal_histname = f'total-deltay'
            if not deltatotal_histname in self.hists:
                title = f'total US y-barycentre - track y-position;#Delta y [cm];Counts'
                self.hists[deltatotal_histname] = ROOT.TH1F(deltatotal_histname, title, 100, -50, 50)
            self.hists[deltatotal_histname].Fill(y_barycentre - yEx)    
            # Fill dy hist per plane
            deltaplane_histname = f'plane{plane}-deltay'
            if not deltaplane_histname in self.hists:
                title = '#splitline{y-barycentre - track y-position}{plane '+str(plane+1)+'};#Delta y [cm];Counts'
                self.hists[deltaplane_histname] = ROOT.TH1F(deltaplane_histname, title, 100, -50, 50)
            self.hists[deltaplane_histname].Fill(y_barycentre - yEx)

            # Fill dy v ypos hist in system
            total_histname = f'total-dyvy'
            if not total_histname in self.hists:
                title = 'total US (y-barycentre - y-position);#Delta(y_{B}, y_{track}) [cm];y_{track} [cm];Counts'
                self.hists[total_histname] = ROOT.TH2F(total_histname, title, 100, -50, 50, 100, -10, 90)
            self.hists[total_histname].Fill(y_barycentre - yEx, yEx)
            # Fill dy v ypos hist per plane
            plane_histname = f'plane{plane}-dyvy'
            if not plane_histname in self.hists:
                title = '#splitline{(y-barycentre - y-position)}{plane '+str(plane+1)+'};#Delta(y_{B}, y_{track}) [cm];y_{track} [cm];Counts'
                self.hists[plane_histname] = ROOT.TH2F(plane_histname, title, 100, -50, 50, 100, -10, 90)
            self.hists[plane_histname].Fill(y_barycentre - yEx, yEx )
            
            # Fill xEx hist
            xEx_histname = f'xEx'
            if not xEx_histname in self.hists:
                title = f'Track x-position;x-position [cm];Counts'
                self.hists[xEx_histname] = ROOT.TH1F(xEx_histname, title, 120, -110, 10)
            self.hists[xEx_histname].Fill(xEx)    
            # Fill yEx hist
            yEx_histname = f'yEx'
            if not yEx_histname in self.hists:
                title = f'Track y-position;y-position [cm];Counts'
                self.hists[yEx_histname] = ROOT.TH1F(yEx_histname, title, 100, -10, 90)
            self.hists[yEx_histname].Fill(yEx) 



    def ScifiCorrectedTimes(self, hit):
        detID=hit.GetDetectorID()

        correctedtimes=self.muAna.GetCorrectedTimes(hit, self.tw.Ex.x(), mode='aligned')
    
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

            self.hists[SiPMtime].Fill(time)

    def XTHists(self, hit):

        correctedtimes, qdcs = self.muAna.GetCorrectedTimes(hit, self.tw.Ex.x(), mode='aligned'), hit.GetAllSignals()
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

    def WriteOutReconstructionHistograms(self):
        if not self.simulation: d = f'{self.outpath}splitfiles/run{self.runNr}/SystemAlignment/'
        else: d = f'{self.outpath}SystemAlignment/'
        os.makedirs(d, exist_ok=True)
        outfilename=d+f'muonreconstruction_{self.options.nStart}.root'

        if os.path.exists(outfilename): outfile=ROOT.TFile.Open(outfilename, 'recreate')
        else: outfile=ROOT.TFile.Open(outfilename, 'create')

        for hname in self.hists:
            
            hist=self.hists[hname]
            if hname in ('yEx', 'xEx', 'reft'):
                outfile.WriteObject(hist, hname, 'kOverwrite')
                continue

            x, proj = hname.split('-')

            if not hasattr(outfile, proj): folder=outfile.mkdir(proj)
            else: folder=outfile.Get(proj)
            
            folder.cd()
            hist.Write(hname, 2) # The 2 means it will overwrite a hist of the same name            

        outfile.Close()
        print(f'{len(self.hists)} histograms saved to {outfilename}')        

    def WriteOutHistograms(self):

        if not self.simulation: d = f'{self.outpath}splitfiles/run{self.runNr}/SystemAlignment/'
        else: d = f'{self.outpath}SystemAlignment/'


        os.makedirs(d, exist_ok=True)
        outfilename=d+f'SystemAlignment_{self.options.nStart}.root'
        
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
