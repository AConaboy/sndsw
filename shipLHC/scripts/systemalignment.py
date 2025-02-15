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
        self.muontrackIDs = [0,1]

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

        self.freq=160.316E6
        self.TDC2ns=1E9/self.freq

        self.sides=('left', 'right')

        if options.XT:
            from itertools import combinations
            self.vetocombinations={'left':list(combinations(range(8), 2)), 'right':list(combinations(range(8, 16), 2))}
            self.UScombinations={'left':list(combinations([0,1,3,4,6,7], 2)), 'right':list(combinations([8,9,11,12,14,15], 2))}

            self.SiPMcombinations = {1:{'left':list(combinations(range(8), 2)), 'right':list(combinations(range(8, 16), 2))}, 
            2:{'left':list(combinations([0,1,3,4,6,7], 2)), 'right':list(combinations([8,9,11,12,14,15], 2))}
            }

        if self.tw.mode=='reconstructmuonposition':
            from buildhistograms import MakeReconstructMuonHistograms
            self.hists = MakeReconstructMuonHistograms()

        else: self.hists=tw.hists

    def FillSiPMHists(self, hit):
        detID=hit.GetDetectorID()
        s,p,b=self.muAna.parseDetID(detID)

        if not self.simulation: 
            correctedtimes=self.muAna.GetCorrectedTimes(hit, MuFilter=self.tw.MuFilter,x=self.tw.Ex[p].x(), mode='aligned')
        else: # simulation
            correctedtimes=hit.GetAllTimes()

        for ch in correctedtimes:
            SiPM, time = ch
            d = self.muAna.alignmentparameters[f'{detID}_{SiPM}']
            # print(f'{detID}_{SiPM}: time = {time}, d = {d[0]}')
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
        nleft, nright = self.muAna.GetnFiredSiPMs(hit)

        s, p, b = self.muAna.parseDetID(detID)
        self.tw.MuFilter.GetPosition(detID, self.A, self.B)
        x_midpoint = 0.5 * (self.A.x() + self.B.x())        

        if not self.simulation:
            correctedtimes = self.muAna.GetCorrectedTimes(hit, MuFilter=self.tw.MuFilter,x=0, mode='aligned')
            tofcorrectedtimes = self.muAna.GetCorrectedTimes(hit, MuFilter=self.tw.MuFilter,x=self.tw.Ex[p].x(), mode='aligned')
        else: 
            correctedtimes = hit.GetAllTimes()
            tofcorrectedtimes = self.muAna.GetCorrectedTimes(hit, x=self.tw.Ex[p].x(), mode='tof')

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
        averagetime = 1/2 * sum(averages.values())
        
        # For determining bar side-wise time resolution
        tofa_averages = {side:(sum(tofaligned_times[side].values()) / len(tofaligned_times[side])) for side in tofaligned_times}
        
        averagebartimehistname=f'averagetime_{detID}_aligned'
        if not averagebartimehistname in self.hists:
            title=self.subsystemdict[s]+' plane '+str(p+1)+' bar '+str(b+1)+' average of aligned times from each side as a function of x_{predicted};x_{predicted} [cm];'+'#frac{1}{2}#times(t^{DS,TW}_{left}+t^{DS,TW}_{right}) [ns];Counts'
            histformat = self.tw.histformatting['averagetime']['aligned']
            self.hists[averagebartimehistname]=ROOT.TH2F(averagebartimehistname, title, *histformat[0], *histformat[1])
        self.hists[averagebartimehistname].Fill(self.tw.Ex[p].x(), averagetime)

        for side in ('left', 'right'):
            averagebarsidetimehistname=f'sidetime_{detID}-{side}_aligned'
            SiPMcut_averagebarsidetimehistname=f'sidetime_{detID}-{side}_aligned-SiPMcut'
            if not averagebarsidetimehistname in self.hists:
                title=self.subsystemdict[s]+' plane '+str(p+1)+' bar '+str(b+1)+' average aligned + ToF corr time from '+side+' side as a function of x_{predicted};x_{predicted} [cm];'+ side+' side average of t^{DS,TW}_{'+side+'} [ns];Counts'
                histformat=self.tw.histformatting['sidetime']['aligned']
                self.hists[averagebarsidetimehistname]=ROOT.TH2F(averagebarsidetimehistname, title, *histformat[0], *histformat[1])
            self.hists[averagebarsidetimehistname].Fill(self.tw.Ex[p].x(), tofa_averages[side])
            
            if not SiPMcut_averagebarsidetimehistname in self.hists:
                title=self.subsystemdict[s]+' plane '+str(p+1)+' bar '+str(b+1)+' average aligned + ToF corr time from '+side+' side as a function of x_{predicted};x_{predicted} [cm];'+ side+' side average of t^{DS,TW}_{'+side+'} [ns];Counts'
                histformat = self.tw.histformatting['sidetime']['aligned']
                self.hists[SiPMcut_averagebarsidetimehistname]=ROOT.TH2F(SiPMcut_averagebarsidetimehistname, title, *histformat[0], *histformat[1])
            if self.muAna.AllLiveSiPMs(hit):
                self.hists[SiPMcut_averagebarsidetimehistname].Fill(self.tw.Ex[p].x(), tofa_averages[side])

    def ReconstructMuonPosition(self, hits):

<<<<<<< HEAD
        barycentres=self.muAna.GetBarycentres(hits)
>>>>>>> 6209bb093 (Updates to fix chi2 in data)
=======
        barycentres=self.muAna.GetBarycentres(hits, MuFilter=self.tw.MuFilter)
>>>>>>> 6fe2275fe (no idea what I'm doing)
        x_methods_dict={mode:self.muAna.GetOverallXBarycentre(barycentres, mode=mode) for mode in ('relQDC', 'maxQDC')}

        for plane in barycentres:
            if len(barycentres[plane])==0: continue
            x_data, y_data = barycentres[plane].values()

            if self.tw.hasTrack: 
                
                xEx, yEx, zEx = self.tw.muAna.GetExtrapolatedPosition(plane)

            # y-histograms
            self.hists[f'lambda_y-plane{plane}'].Fill(y_data['lambda_y'])

            if self.tw.hasTrack: 
                dy = y_data['yB'] - yEx

                self.hists[f'dyB-plane{plane}'].Fill(dy)
                self.hists[f'dyB-total'].Fill(dy)

                self.hists[f'dyBvEy-plane{plane}'].Fill(dy, yEx)
                self.hists[f'dyBvEy-total'].Fill(dy, yEx)

                self.hists[f'lambda_yvEy-plane{plane}'].Fill(y_data["lambda_y"], yEx)

            for key,x_method in x_methods_dict.items():
                        
                x_val=x_method[plane]['lambda_x']
                self.hists[f'lambda_x_{key}-plane{plane}'].Fill(x_val)

            if not self.tw.hasTrack: return

            for key,x_method in x_methods_dict.items():                
                
                for x in ["dxL", "dxR", "dxB"]:
                    
                    if x in ('dxL', 'dxR'): x_val = x_method[plane][x][0]-xEx
                    else: x_val = x_method[plane][x]-xEx

                    self.hists[f'{x}_{key}-plane{plane}'].Fill(x_val)
                    self.hists[f'{x}_{key}-total'].Fill(x_val)

                    self.hists[f'{x}_{key}vEx-plane{plane}'].Fill(x_val, xEx)
                    self.hists[f'{x}_{key}vEx-total'].Fill(x_val, xEx)

                    self.hists[f'lambda_x_{key}vEx-plane{plane}'].Fill(x_val, xEx)

            for detID in x_data:
                relQDC=x_data[detID]['relQDC']
                self.hists[f'relQDC-plane{plane}'].Fill(relQDC)
                self.hists[f'relQDC-total'].Fill(relQDC)

    def ScifiCorrectedTimes(self, hit):
        detID=hit.GetDetectorID()

        correctedtimes=self.muAna.GetCorrectedTimes(hit, self.tw.Ex[p].x(), mode='aligned')
    
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
                histformat=self.tw.histformatting['dtvxpred']['corrected']
                self.hists[SiPMtime]=ROOT.TH1F(SiPMtime,fulltitle, *histformat[1])

            self.hists[SiPMtime].Fill(time)

    def XTHists(self, hit):

        detID = hit.GetDetectorID()
        s,p,b = self.muAna.parseDetID(detID)

        correctedtimes = self.muAna.GetCorrectedTimes(hit, MuFilter=self.tw.MuFilter,x=self.tw.Ex[p].x(), mode='aligned')
        qdcs = hit.GetAllSignals()
        times = {'left':{}, 'right':{}}
        
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
            histformat=self.tw.histformatting['timingxt']['aligned']
            if self.timealignment=='old': self.hists[xrefthistname]=ROOT.TH2F(xrefthistname,fulltitle, *histformat[0], *histformat[1])
            else: self.hists[xrefthistname]=ROOT.TH2F(xrefthistname,fulltitle,*histformat[0], *histformat[1])

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
                    histformat=self.tw.histformatting['timingxt']['aligned']
                    if self.timealignment=='old': self.hists[xthistname]=ROOT.TH2F(xthistname,fulltitle, *histformat[0], *histformat[1])
                    else: self.hists[xthistname]=ROOT.TH2F(xthistname,fulltitle, *histformat[0], *histformat[1])

                self.hists[xthistname].Fill(time_i, time_j)
                self.hists[xrefthistname].Fill(time_i, self.tw.reft)

    def MultiMuonReconstruction(self, hits):

        samebar, leftID, rightID = self.GetMuonTracks()

        if not samebar:
            print(f'Muons not in same bar, skipping event {self.tw.M.EventNumber}')
            return 

        barycentres=self.muAna.GetBarycentres(hits, MuFilter=self.tw.MuFilter)
        x_barycentres = self.muAna.GetOverallXBarycentre(barycentres, mode='maxQDC')

        left_muon = self.tw.M.eventTree.MCTrack[leftID]
        right_muon = self.tw.M.eventTree.MCTrack[rightID]

        for plane in barycentres:
            if len(barycentres[plane])==0: continue
            
            # Get detID in plane where both muons are
            muons_detID = self.GetMultiMuonBar(plane)
            if muons_detID== False: continue

            x_data, y_data = barycentres[plane].values()

            # Get the track extrapolations for the left muon and right muon
            lefttrack_pos = self.ExtrapolateMCTrack(leftID, plane)
            righttrack_pos = self.ExtrapolateMCTrack(rightID, plane)
            
            xL = x_data[muons_detID]['xL']
            xR = x_data[muons_detID]['xR']
            
            # Residuals between measurements of xL/xR and the muon positions
            dxL = lefttrack_pos[0] - xL 
            dxR = righttrack_pos[0] - xR 

            dxLvtrackx_histname = f'dxL-plane{plane}'
            if not dxLvtrackx_histname in self.hists:
                title = 'Residual between x_{left} and left muon track position for plane '+str(plane+1)+';#Delta(x_{left}, left muon track x) [cm];left muon track x [cm]'
                self.hists[dxLvtrackx_histname] = ROOT.TH2F(dxLvtrackx_histname, title, 50, -25, 25, 90, -80, 10)
            self.hists[dxLvtrackx_histname].Fill(dxL, lefttrack_pos[0])

            if not dxRvtrackx_histname in self.hists:
                title = 'Residual between x_{right} and right muon track position for plane '+str(plane+1)+';#Delta(x_{right}, right muon track x) [cm];right muon track x [cm]'
                self.hists[dxRvtrackx_histname] = ROOT.TH2F(dxRvtrackx_histname, title, 50, -25, 25, 90, -80, 10)                
            self.hists[dxRvtrackx_histname].Fill(dxR, righttrack_pos[0])

    def ExtrapolateMCTrack(self, trackID, plane):
        track = self.tw.M.eventTree.MCTrack[trackID]

        self.tw.MuFilter.GetPosition(int(f'2{plane}000'), self.tw.A, self.tw.B)
        planez = 0.5*(self.tw.A.z() + self.tw.B.z())

        startx, starty, startz = track.GetStartX(), track.GetStartY(), track.GetStartZ()
        px, py, pz = track.GetPx(), track.GetPy(), track.GetPz()

        pos = [startx+px/pz*(planez-startz), starty+py/pz*(planez-startz), planez]
        return pos

    def GetMultiMuonBar(self, plane):

        tracks = {i:self.tw.M.eventTree.MCTrack[i] for i in [0,1]}
        
        start_xs = {i:tracks[i].GetStartX() for i in self.muontrackIDs}
        start_ys = {i:tracks[i].GetStartY() for i in self.muontrackIDs}
        start_zs = {i:tracks[i].GetStartX() for i in self.muontrackIDs}
        
        pxs = {i:tracks[i].GetPx() for i in self.muontrackIDs}
        pys = {i:tracks[i].GetPy() for i in self.muontrackIDs}
        pzs = {i:tracks[i].GetPz() for i in self.muontrackIDs}
        
        self.tw.MuFilter.GetPosition(int(f'2{plane}000'), self.tw.A, self.tw.B)
        plane_z = 0.5*(self.tw.A.z() + self.tw.B.z())

        x_planez = {i:start_xs[i]+pxs[i]/pzs[i]*(plane_z - start_zs[i]) for i in self.muontrackIDs}
        y_planez = {i:start_ys[i]+pys[i]/pzs[i]*(plane_z - start_zs[i]) for i in self.muontrackIDs}

        detIDs = {i:self.tw.nav.FindNode(x_planez[i], x_planez[i], plane_z).GetName() for i in self.muontrackIDs}
        if detIDs[1] != detIDs[0]: return False 
        elif detIDs[1].find('volMuUpstreamBar'): 
            print(f'Not in HCAL bar: {detIDs[1]}')
            return False
        else: 
            print(detIDs[1])
            return int(detIDs[1].split('_')[-1])

    def GetMuonTracks(self):

        tracks = {i:self.tw.M.eventTree.MCTrack[i] for i in self.muontrackIDs}
        
        start_xs = {i:tracks[i].GetStartX() for i in self.muontrackIDs}
        start_ys = {i:tracks[i].GetStartY() for i in self.muontrackIDs}
        start_zs = {i:tracks[i].GetStartX() for i in self.muontrackIDs}
        
        pxs = {i:tracks[i].GetPx() for i in self.muontrackIDs}
        pys = {i:tracks[i].GetPy() for i in self.muontrackIDs}
        pzs = {i:tracks[i].GetPz() for i in self.muontrackIDs}
        
        # Extrapolate tracks to HCAL plane 1
        self.tw.MuFilter.GetPosition(20000, self.tw.A, self.tw.B)
        HCAL1z = 0.5*(self.tw.A.z() + self.tw.B.z())
        # Below is the x pos of each track at HCAL1-z
        x_HCAL1z = {i:start_xs[i]+pxs[i]/pzs[i]*(HCAL1z - start_zs[i]) for i in self.muontrackIDs}
        y_HCAL1z = {i:start_ys[i]+pys[i]/pzs[i]*(HCAL1z - start_zs[i]) for i in self.muontrackIDs}

        # Check if muons in the same bar
        if self.tw.nav.FindNode(x_HCAL1z[0], y_HCAL1z[0], HCAL1z).GetName() != self.tw.nav.FindNode(x_HCAL1z[1], y_HCAL1z[1], HCAL1z).GetName(): samebar=False 
        else: samebar=True

        rightID = min(x_HCAL1z, key=x_HCAL1z.get)
        leftID = 1 if rightID==2 else 2
        return samebar, leftID, rightID

    def WriteOutReconstructionHistograms(self):
        if not self.simulation: 
            d = f'{self.outpath}splitfiles/run{self.runNr}/{self.tw.mode}/'
            os.makedirs(d, exist_ok=True)
            
            outfilename=d+f'muonreconstruction_{self.options.nStart}.root'

        elif self.simulation and self.options.simMode=='muonDIS': 
            d = f'{self.outpath}{self.tw.mode}/'
            os.makedirs(d, exist_ok=True)

            preamble, key1, key2, end = self.options.fname.split('_')
            outfilename=d+f'muonreconstruction_{key1}_{key2}.root'

        elif self.simulation and self.options.simMode=='neutrino': 
            d = f'{self.outpath}{self.tw.mode}/'
            os.makedirs(d, exist_ok=True)            

            simulation_details, filename = self.options.fname.split('/')
            key=filename.split('_')[1]
            outfilename=d+f'muonreconstruction_{key}.root'

        elif self.simulation and self.options.simMode == 'neutralhadron':
            d = f'{self.outpath}{self.tw.mode}/'
            os.makedirs(d, exist_ok=True)

            particle_type, Emin, Emax, key = self.options.fname.split('_')[3:7]
            outfilename=d+f'muonreconstruction_{particle_type}_{Emin}_{Emax}_{key}.root'

        elif self.simulation and self.options.simMode == 'passingmuon':
            d = f'{self.outpath}{self.tw.mode}/'
            os.makedirs(d, exist_ok=True)

            keys=self.options.fname.split('_')[1:3]
            outfilename=d+f'muonreconstruction_{keys[0]}_{keys[1]}.root'

        elif self.simulation and self.options.simMode == 'nue':
            print(f'Not implemented nue saving protocol! ')
      
        if os.path.exists(outfilename): outfile=ROOT.TFile.Open(outfilename, 'recreate')
        else: outfile=ROOT.TFile.Open(outfilename, 'create')

        for hname in self.hists:
            
            hist=self.hists[hname]
            if hname in ('yEx', 'xEx', 'reft'):
                outfile.WriteObject(hist, hname, 'kOverwrite')
                continue

            key, plane = hname.split('-')

            if not hasattr(outfile, key): folder=outfile.mkdir(key)
            else: folder=outfile.Get(key)
            
            folder.cd()
            hist.Write(hname, 2) # The 2 means it will overwrite a hist of the same name            

        outfile.Close()
        print(f'{len(self.hists)} histograms saved to {outfilename}')        

    def WriteOutHistograms(self):

        if not self.simulation: 
            d = f'{self.outpath}splitfiles/run{self.runNr}/{self.tw.mode}/'
            os.makedirs(d, exist_ok=True)
            outfilename=d+f'systemalignment_{self.options.nStart}.root'

        elif self.simulation and self.options.simMode=='muonDIS': 
            d = f'{self.outpath}{self.tw.mode}/'
            os.makedirs(d, exist_ok=True)

            preamble, key1, key2, end = self.options.fname.split('_')
            outfilename=d+f'systemalignment_{key1}_{key2}.root'

        elif self.simulation and self.options.simMode=='neutrino': 
            d = f'{self.outpath}{self.tw.mode}/'
            os.makedirs(d, exist_ok=True)            

            simulation_details, filename = self.options.fname.split('/')
            key=filename.split('_')[1]
            outfilename=d+f'systemalignment_{key}.root'

        elif self.simulation and self.options.simMode == 'neutralhadron':
            d = f'{self.outpath}{self.tw.mode}/'
            os.makedirs(d, exist_ok=True)

            particle_type, Emin, Emax, key = self.options.fname.split('_')[3:7]
            outfilename=d+f'systemalignment_{particle_type}_{Emin}_{Emax}_{key}.root'

        elif self.simulation and self.options.simMode == 'passingmuon':
            keys=self.options.fname.split('_')[1:3]
            outfilename=d+f'systemalignment_{keys[0]}_{keys[1]}.root'

        elif self.simulation and self.options.simMode == 'nue':
            print(f'Not implemented nue saving protocol! ')        
        
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
