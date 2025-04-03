#!/usr/bin/env python
import ROOT, csv, os, json
import numpy as np
from pathlib import Path
from HCALTools import HCALTools 
ROOT.gInterpreter.ProcessLine('#include "/afs/cern.ch/user/a/aconsnd/sndsw/analysis/tools/sndSciFiTools.h"')

class ShowerProfiles(object):

    def __init__(self, options, tw):
       
        self.tw=tw
        self.simulation=tw.simulation
        self.runNr = tw.runNr
        self.muAna = tw.muAna
        self.timealignment=tw.timealignment
        self.TWCorrectionRun=tw.TWCorrectionRun

        ### If no time-walk correction run is provided. Set the default correction run depending on time alignment of the data set

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

        self.MuFilter = tw.MuFilter
        self.barlengths = self.muAna.BuildBarLengths(self.MuFilter)

        self.freq=160.316E6
        self.TDC2ns=1E9/self.freq

        self.sides=('left', 'right')

        from buildhistograms import MakeShowerProfilesHistograms
        self.hists = MakeShowerProfilesHistograms()
        
        self.A, self.B = ROOT.TVector3(), ROOT.TVector3()
        
        if options.signalpartitions: self.Loadnumuevents()
        
        if options.numuStudy or options.nueStudy: self.nuStudy=True 
        else: self.nuStudy=False 

        if self.nuStudy: self.all_barycentres={}
    
        self.barycentres={}
        self.clusters={}
        self.scifidata={}
        self.data={}
        
        self.hcalTools = HCALTools(self.muAna, self.tw.MuFilter)
        self.InstanceCuts()

        if self.simulation:

            self.column_names=['EventNumber', 'scaled_weight',
            'hasMuon','interactionWall',
            'lambdax0','lambdax1', 'lambdax2','lambdax3', 'lambdax4', 
            'lambday0','lambday1', 'lambday2','lambday3', 'lambday4',
            'relQDC0','relQDC1','relQDC2','relQDC3','relQDC4',
            'planeQDC0','planeQDC1','planeQDC2','planeQDC3','planeQDC4',
            'hasTrack'
            ]
            [self.column_names.append(i) for i in self.cuts]

            # Instance list of lists to store data
            # Estimated to be about 32 MB/500k events using float32
            self.eventdata=[self.column_names]
        
    def InstanceCuts(self):
    
        listofcuts=['veto', 'avgSFchannel', 'avgDSchannel', 'SF1', 'SF2', '2consecSFplanes', 'allHCALplanes']
        
        self.cuts = {"vetoCut": ROOT.snd.analysis_cuts.vetoCut(self.tw.M.eventTree),
            "avgSFchannel":ROOT.snd.analysis_cuts.avgSciFiFiducialCut(200, 1200, 300, 128*12-200, self.tw.M.eventTree),
            "avgDSchannel":ROOT.snd.analysis_cuts.avgDSFiducialCut(70, 105, 10, 50, self.tw.M.eventTree),
            "SF1":ROOT.snd.analysis_cuts.sciFiStationCut(0., ROOT.std.vector('int')([1]), self.tw.M.eventTree),
            "SF2":ROOT.snd.analysis_cuts.sciFiStationCut(0., ROOT.std.vector('int')([2]), self.tw.M.eventTree),
            "consec2SFplanes":ROOT.snd.analysis_cuts.minSciFiConsecutivePlanes(self.tw.M.eventTree),
            "allHCALplanes":ROOT.snd.analysis_cuts.DSActivityCut(self.tw.M.eventTree)
                }

    def Loadnumuevents(self):
        numusignalevent_filepath = '/afs/cern.ch/work/a/aconsnd/numusignalevents.csv'
        with open(numusignalevent_filepath, 'r') as f:
            reader=csv.reader(f)
            nu_mu_data=[r for r in reader]
        self.nu_mu_events={int(x[0]):(int(x[1]), int(x[2])) for x in nu_mu_data[1:]}        

    def FillHists(self, hits):
        
        """ 
        Hits passed here will be in events with a DS track
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
            
            frac_n = ( len(found_large_SiPMs) + len(found_small_SiPMs) ) / 16
            self.hists[fractionSiPMs_histname].Fill(frac_n) 

            if self.tw.options.SmallSiPMcheck:
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

    def ExtractScifiData(self, hits):
        self.runId = self.tw.M.eventTree.EventHeader.GetRunId()
        self.scifidata[self.runId] = {}
        for idx, hit in enumerate(hits):
            detID = hit.GetDetectorID()
            
            self.tw.M.Scifi.GetSiPMPosition(detID, self.A,self.B)
            if hit.isVertical(): pos = 0.5*(self.A[0] + self.B[0])
            else: pos = 0.5*(self.A[1] + self.B[1])

            signal, ctime = hit.GetSignal(), self.tw.M.Scifi.GetCorrectedTime(detID, hit.GetTime(), 0)
            self.scifidata[self.runId][idx]=[detID, pos, signal, ctime]

    def ScifiClusterInfo(self, clusters):
        self.runId = self.tw.M.eventTree.EventHeader.GetRunId()
        self.clusters[self.runId] = {}
        interactionWall = self.GetInteractionWall(self.runId)
        seedpos = self.GetSeedPos(self.runId)

        for idx, cl in enumerate(clusters):
            first=cl.GetFirst()
            cl.GetPosition(self.A, self.B)
            avg_x, avg_y = 0.5*(self.A.x()+self.B.x()), 0.5*(self.A.y()+self.B.y())
            N, totalQDC = cl.GetN(), cl.GetEnergy()
            x0,y0 = seedpos
            self.clusters[self.runId][idx]=[interactionWall,first, avg_x, avg_y, N, totalQDC, x0, y0]

    def ShowerDirection(self, hits, scifi_hits):
        
        """
        validity of hits is ensured in muAna.GetPlaneData(hits, mufilter) ! 
        """

        self.EventNumber = self.tw.M.EventNumber
        self.interactionWall = self.hcalTools.GetInteractionWall(scifi_hits)
        self.barycentres = self.muAna.GetBarycentres(hits, MuFilter=self.tw.MuFilter)
        self.x_methods_dict = {mode:self.muAna.GetOverallXBarycentre(self.barycentres, mode=mode) for mode in ('relQDC', 'maxQDC')}        
        
        # For numu/nue study, 
        if self.nuStudy: 
            self.runId = self.tw.M.eventTree.EventHeader.GetRunId()
            if int(self.runId) in self.all_barycentres: 
                print(f'run number {self.runId} already in all barycentres dict')
            self.all_barycentres[self.runId] = self.barycentres 

        """
        1. Plot shower barycentre in y (x) by using a QDC weighting, aligned timing
        2. First consider all hits 
        3. For each event, plot the barycentre in x and y at each plane
        """

        self.FillBarycentrePlots(hits)

        # Write out all SiPM information when looking at the numu data
        if self.nuStudy: 
            self.RecordAllSiPMdata(hits)

        self.pass_cuts = {k:v.passCut() for k,v in self.cuts.items()}  # should work with sim or not
        if self.simulation: self.StoreSimulationDetails()

    def FillBarycentrePlots(self, hits):

        self.lambda_x_dict = {i:np.nan for i in range(5)}
        self.lambda_y_dict = {i:np.nan for i in range(5)}
        self.totalQDC_dict = {i:np.nan for i in range(5)}
        self.relQDC_dict = {i:np.nan for i in range(5)}

        for plane in self.barycentres:
            # All planes are instanced as empty dictionaries
            # to retain the correct order for passing to the BDT
            if len(self.barycentres[plane])==0: continue 

            x_data, y_data = self.barycentres[plane].values()

            # Get bar max QDC and relQDC of that bar for this plane
            maxQDC_detID = max(x_data, key=lambda detID: x_data[detID]['barQDC'])
            maxQDC = x_data[maxQDC_detID]['barQDC']
            totalQDC = sum(x_data[detID]['barQDC'] for detID in x_data)
            relQDC = maxQDC / totalQDC

            # print(f'relQDC:{relQDC}, totalQDC:{totalQDC}')

            if self.tw.hasTrack: 
                xEx, yEx, zEx = self.tw.muAna.GetExtrapolatedPosition(plane)

            # y-histograms
            if self.simulation: 
                self.hists[f'lambda_y-plane{plane}'].Fill(y_data['lambda_y'], self.tw.scaled_event_weight)
            else: 
                self.hists[f'lambda_y-plane{plane}'].Fill(y_data['lambda_y'])
            self.lambda_y_dict[plane]=y_data['lambda_y']

            if self.tw.hasTrack: 
                xEx, yEx, zEx = self.tw.muAna.GetExtrapolatedPosition(plane)
                dy = y_data['yB'] - yEx

                if self.simulation:
                    # Make y_b histograms
                    self.hists[f'dyB-plane{plane}'].Fill(dy, self.tw.scaled_event_weight)
                    self.hists[f'dyB-total'].Fill(dy, self.tw.scaled_event_weight)

                    self.hists[f'dyBvEy-plane{plane}'].Fill(dy, yEx, self.tw.scaled_event_weight)
                    self.hists[f'dyBvEy-total'].Fill(dy, yEx, self.tw.scaled_event_weight)

                    self.hists[f'lambda_yvEy-plane{plane}'].Fill(y_data["lambda_y"], yEx, self.tw.scaled_event_weight)
                else: 
                    # Make y_b histograms
                    self.hists[f'dyB-plane{plane}'].Fill(dy)
                    self.hists[f'dyB-total'].Fill(dy)

                    self.hists[f'dyBvEy-plane{plane}'].Fill(dy, yEx)
                    self.hists[f'dyBvEy-total'].Fill(dy, yEx)                    
                    self.hists[f'lambda_yvEy-plane{plane}'].Fill(y_data["lambda_y"], yEx)

            for key,x_method in self.x_methods_dict.items():
                        
                x_val=x_method[plane]['lambda_x']
                
                if self.simulation: 
                    self.hists[f'lambda_x_{key}-plane{plane}'].Fill(x_val, self.tw.scaled_event_weight)
                    self.hists[f'lambda_x_{key}vQDC-plane{plane}'].Fill(x_val, maxQDC, self.tw.scaled_event_weight)
                else: 
                    self.hists[f'lambda_x_{key}-plane{plane}'].Fill(x_val)
                    self.hists[f'lambda_x_{key}vQDC-plane{plane}'].Fill(x_val, maxQDC)

            self.lambda_x_dict[plane] = self.x_methods_dict['maxQDC'][plane]['lambda_x']
            self.totalQDC_dict[plane] = totalQDC
            self.relQDC_dict[plane] = relQDC 

            if not self.tw.hasTrack: return

            for key,x_method in self.x_methods_dict.items():                
                
                for x in ["dxL", "dxR", "dxB"]:
                    
                    if x in ('dxL', 'dxR'): x_val = x_method[plane][x][0]-xEx
                    else: x_val = x_method[plane][x]-xEx

                    if self.simulation:
                        self.hists[f'{x}_{key}-plane{plane}'].Fill(x_val, self.tw.scaled_event_weight)
                        self.hists[f'{x}_{key}-total'].Fill(x_val, self.tw.scaled_event_weight)

                        self.hists[f'{x}_{key}vEx-plane{plane}'].Fill(x_val, xEx, self.tw.scaled_event_weight)
                        self.hists[f'{x}_{key}vEx-total'].Fill(x_val, xEx, self.tw.scaled_event_weight)

                        self.hists[f'lambda_x_{key}vEx-plane{plane}'].Fill(x_val, xEx, self.tw.scaled_event_weight)
                    else:
                        self.hists[f'{x}_{key}-plane{plane}'].Fill(x_val)
                        self.hists[f'{x}_{key}-total'].Fill(x_val)

                        self.hists[f'{x}_{key}vEx-plane{plane}'].Fill(x_val, xEx)
                        self.hists[f'{x}_{key}vEx-total'].Fill(x_val, xEx)

                        self.hists[f'lambda_x_{key}vEx-plane{plane}'].Fill(x_val, xEx) 

            for detID in x_data:
                # relQDC=x_data[detID]['relQDC']
                if self.simulation:
                    self.hists[f'relQDC-plane{plane}'].Fill(relQDC, self.tw.scaled_event_weight)
                    self.hists[f'relQDC-total'].Fill(relQDC, self.tw.scaled_event_weight)
                else:
                    self.hists[f'relQDC-plane{plane}'].Fill(relQDC)
                    self.hists[f'relQDC-total'].Fill(relQDC)

    def StoreSimulationDetails(self):

        output_dict = {k:None for k in self.column_names}

        if self.tw.options.simMode=='neutrino':
            output_dict['flav'] = self.GetNeutrinoIntType(self.tw.M.eventTree)

        eventHasMuon=self.hcalTools.OutgoingMuon(self.tw.M.eventTree, self.tw.options.simMode)
        output_dict['EventNumber'] = self.EventNumber
        output_dict['scaled_weight'] = self.tw.scaled_event_weight
        output_dict['hasMuon'] = eventHasMuon # truth of whether there is a final state muon
        output_dict['hasTrack'] = self.tw.hasTrack # is track found? 
        output_dict['interactionWall'] = self.interactionWall

        # print(self.relQDC_dict, self.totalQDC_dict)

        for i in range(5):
            output_dict[f'lambdax{i}'] = np.float32(self.lambda_x_dict[i])
            output_dict[f'lambday{i}'] = np.float32(self.lambda_y_dict[i])
            output_dict[f'relQDC{i}'] = np.float32(self.relQDC_dict[i])
            output_dict[f'planeQDC{i}'] = np.float32(self.totalQDC_dict[i])

        for cutname in self.pass_cuts: output_dict[cutname] = self.pass_cuts[cutname]

        output = [output_dict[key] for key in self.column_names]
        
        # Not writing out events where all lambda values are nan
        all_lambdaxs_nan = all([np.isnan(output_dict[f'lambdax{i}']) for i in range(5) ])
        all_lambdays_nan = all([np.isnan(output_dict[f'lambday{i}']) for i in range(5) ])
        if all_lambdaxs_nan and all_lambdays_nan: return
        
        self.eventdata.append(output)
        
    def WriteOutSimulationDetails(self):
        d, outfilename_noExt = self.GetOutFileName()
        datafilename = d+outfilename_noExt+'.csv'
        
        with open(datafilename, 'w', newline='') as f:
            writer=csv.writer(f)
            # writer.writerow(self.column_names)  # Write the header
            writer.writerows(self.eventdata)
        print(f'Simulation details written to {datafilename}')

    def GetBarycentre(self, plane, proj):

        if proj=='x': 
            if plane not in self.xbarycentres: return np.nan
            if not 'dxB' in self.xbarycentres[plane]: return np.nan
            b=self.xbarycentres[plane]['dxB']

        elif proj=='y': 
            if plane not in self.barycentres: return np.nan
            if not 'y-barycentre' in self.barycentres[plane]: return np.nan
            if not 'yB' in self.barycentres[plane]['y-barycentre']: return np.nan
            b=self.barycentres[plane]['y-barycentre']['yB']
        return b

    def RecordAllSiPMdata(self, hits):
        
        for hit in hits:
            data={}

            detID = hit.GetDetectorID()
            s,p,b = self.muAna.parseDetID(detID)
            if s==3: continue
            interactionWall = self.GetInteractionWall(self.runId)

            xEx, yEx, zEx = self.muAna.GetExtrapolatedPosition(p)

            qdcs = [(SiPM, qdc) for SiPM,qdc in hit.GetAllSignals()]
            rawtimes = [(SiPM, cc*self.TDC2ns) for SiPM,cc in hit.GetAllTimes()]
            twctimes = self.muAna.GetCorrectedTimes(hit, MuFilter=self.MuFilter, mode='tw')
            atimes = self.muAna.GetCorrectedTimes(hit, MuFilter=self.MuFilter, mode='aligned')
            
            if self.muAna.GetExtrapolatedBarDetID(p) == detID: trackrelated = True
            else: trackrelated = False

            for idx, lst in enumerate([qdcs, rawtimes, twctimes, atimes]): 
                for SiPM, t in lst:
                    if not SiPM in data: data[SiPM] = {} 

                    if idx==0:data[SiPM]['qdc']=t 
                    elif idx==1:data[SiPM]['rawtimes']=t 
                    elif idx==2:data[SiPM]['twctimes']=t 
                    elif idx==3:
                        data[SiPM]['atimes']=t 
                        data[SiPM]['reft']=self.tw.reft
                        
                        if f'{detID}_{SiPM}' in self.muAna.alignmentparameters:
                            data[SiPM]['d'] = self.muAna.alignmentparameters[f'{detID}_{SiPM}']
                            data[SiPM]['trackrelated'] = trackrelated
                            data[SiPM]['cscint'] = self.muAna.cscintvalues[f'{detID}_{SiPM}']

                            # Writing out extrapolated positions
                            data[SiPM]['xEx'] = xEx
                            data[SiPM]['yEx'] = yEx
                            data[SiPM]['Interaction wall'] = interactionWall 
            if not self.runId in self.data: self.data[self.runId]={}
            if not self.tw.M.EventNumber in self.data[self.runId]: self.data[self.runId][self.tw.M.EventNumber]={}
            self.data[self.runId][self.tw.M.EventNumber][detID]=data  

    def GetInteractionWall(self, runNr):
        interactionwall = self.muAna.nu_mu_events[runNr][1]
        return interactionwall
    def GetSeedPos(self, runNr):
        seedpos = self.muAna.nu_mu_events[runNr][2:4]
        return seedpos
    def GetNeutrinoIntType(self, event):

        if not hasattr(event, "MCTrack"):
            print(f'No MCTrack branch. Is this real data?')
            return 

        if event.MCTrack[0].GetPdgCode() == event.MCTrack[1].GetPdgCode():
            i_flav = 0 #NC
        elif abs(event.MCTrack[1].GetPdgCode()) == 11: # electron
            i_flav = 1 #nueCC
        elif abs(event.MCTrack[1].GetPdgCode()) == 13: # muon
            i_flav = 2 #numuCC
        elif abs(event.MCTrack[1].GetPdgCode()) == 15: # tau
            is1Mu = False
            for j_track in range(2, len(event.MCTrack)):
                if event.MCTrack[j_track].GetMotherId() == 1 and abs(event.MCTrack[j_track].GetPdgCode()) == 13:
                    is1Mu = True
                    break
            if is1Mu:
                i_flav = 4 #nutauCC1mu
            else:
                i_flav = 3 #nutauCC0mu 
        else: return 
        
        return i_flav        

    def EvaluateShowerSlopes(self):

        if self.tw.options.simMode=='neutrino':
            i_flav = self.GetNeutrinoIntType(self.tw.M.eventTree)

        # Just make a single hist per proj for evaluating slopes
        for p in ['x', 'y']:
            name=f'{p}z-slope'
            if name in self.hists: self.hists[name].Reset()
            else: 
                title = f'{p}z projection view of Scifi medians and HCAL barycentres;z [cm];{p} [cm]'
                self.hists[name] = ROOT.TH1D(name, title, 250, 250, 500)

        # Add scifi median positions to the plots
        scifi_hitpositions = self.GetScifiMedians()
        for station in scifi_hitpositions:
            
            for orientation in scifi_hitpositions[station]:
                if orientation=='hor':vert=0 
                else: vert=1
                
                zpos=self.tw.zPos['Scifi'][10*station+vert]
                
                data = scifi_hitpositions[station][orientation]
                if len(data)==0: continue # suppresses output for trying to average an empty container

                median = np.median(data)
                err = median/np.sqrt(len(data))

                if orientation=='hor':
                    zbin = self.hists[f'yz-slope'].FindBin(zpos)
                    self.hists[f'yz-slope'].SetBinContent(zbin, median)
                    self.hists[f'yz-slope'].SetBinError(zbin, err)

                elif orientation=='vert':
                    zbin = self.hists[f'xz-slope'].FindBin(zpos)
                    self.hists[f'xz-slope'].SetBinContent(zbin, median)
                    self.hists[f'xz-slope'].SetBinError(zbin, err)                    
        
        # Now add HCAL barycentres to the plot
        for plane in self.barycentres:
            if len(self.barycentres[plane])==0: continue 
            
            x_data, y_data = self.barycentres[plane].values() 

            zPos = self.tw.zPos['MuFilter'][20+plane]
            zbin = self.hists['yz-slope'].FindBin(zPos)

            yB = y_data['yB']
            yBerror = 3 # not sure what else to do

            xB = self.x_methods_dict['maxQDC'][plane]['dxB']
            # Calculate uncertainty
            xLerror = self.x_methods_dict['maxQDC'][plane]['dxL'][1]
            xRerror = self.x_methods_dict['maxQDC'][plane]['dxR'][1]
            xBerror = np.sqrt(xLerror**2 + xRerror**2)

            self.hists['yz-slope'].SetBinContent(zbin, yB)
            self.hists['yz-slope'].SetBinError(zbin, yBerror)

            self.hists['xz-slope'].SetBinContent(zbin, xB)
            self.hists['xz-slope'].SetBinError(zbin, xBerror)

        # Now all data is in the histograms I can determine the slope
        angles={'polar':np.nan, 'azimuthal':np.nan}
        gradients={'x':np.nan, 'y':np.nan}

        for p in ['x', 'y']:
            hist = self.hists[f'{p}z-slope']
            if hist.GetEntries()<3: continue # Can't fit data with less than 3 data points

            # W: use chi-square method, Q: quiet fit, S: return resultPtr, N: don't draw or store the function
            rc = hist.Fit('pol1', 'WQSN') 
            res=rc.Get()
            if not res: continue
            
            gradients[p] = res.Parameter(1)
            gradient_error = res.ParError(1) 
        
        r = np.sqrt(1 + gradients['x']**2 + gradients['y']**2)
        angles['polar'] = np.degrees(np.arccos(1/r))
        angles['azimuthal'] = np.degrees(np.arctan2(gradients['y'], gradients['x']))

        # Fill these angles into a plot
        name = 'polar'
        if not name in self.hists:
            title = '#theta measured with both Scifi and HCAL;Polar angle, #theta_{X} [degrees];Counts'
            self.hists[name] = ROOT.TH1D(name, title, 90, -90, 90)
        if not np.isnan(angles['polar']): self.hists[name].Fill(angles['polar'])
        
        name = 'azimuthal'
        if not name in self.hists:
            title = '#phi measured with both Scifi and HCAL;Azimuthal angle, #phi_{X} [degrees];Counts'
            self.hists[name] = ROOT.TH1D(name, title, 100, -100, 100)
        if not np.isnan(angles['azimuthal']):  self.hists[name].Fill(angles['azimuthal'])

        # thetavphi='thetavphi'
        # if not thetavphi in self.hists:
        #     title = '#theta_{X} v #phi_{X} measured with both Scifi and HCAL;#phi [degrees];#theta [degrees];Counts'
        #     self.hists[thetavphi] = ROOT.TH2D(thetavphi, title, 90, -90, 90, 90, -90, 90)
        # if not np.isnan(angles['x']) and not np.isnan(angles['y']):  self.hists[thetavphi].Fill(angles['x'], angles['y'])

        # Make a plot of the difference between X angle and muon angle
        # Fill these angles into a plot

        if self.tw.hasTrack==False: 
            # print(f'No fitted track in event {self.tw.M.EventNumber}\n')
            return

        slopeX, slopeY = self.tw.mom.x()/self.tw.mom.z(), self.tw.mom.y()/self.tw.mom.z()

        r = np.sqrt(1 + slopeX**2 + slopeY**2)
        muon_polar = np.degrees(np.arccos(1/r))
        muon_azimuthal = np.degrees(np.arctan2(slopeY, slopeX))

        dXmutheta = f'dXmu-polar'
        if not dXmutheta in self.hists:
            title = '#splitline{Difference in polar angle between shower and fitted track, #Delta(#theta_{X}, #theta_{#mu})}{measured with both Scifi and HCAL};Difference in polar angle, #Delta(#theta_{X}, #theta_{#mu}) [degrees];Counts'
            self.hists[dXmutheta] = ROOT.TH1D(dXmutheta, title, 90, -90, 90)
        if not np.isnan(angles['polar']):  self.hists[dXmutheta].Fill(angles['polar'] - muon_polar)
       
        dXmuphi = f'dXmu-azimuthal'
        if not dXmuphi in self.hists:
            title = '#splitline{Difference in azimuthal angle between shower and fitted track, #Delta(#phi_{X}, #phi_{#mu})}{measured with both Scifi and HCAL};Difference in azimuthal angle, #Delta(#phi_{X}, #phi_{#mu}) [degrees];Counts'
            self.hists[dXmuphi] = ROOT.TH1D(dXmuphi, title, 100, -100, 100)
        if not np.isnan(angles['azimuthal']): self.hists[dXmuphi].Fill(angles['azimuthal'] - muon_azimuthal)

        name='scatteringangle'
        if not name in self.hists:
            title = '#splitline{3D scattering angle between shower and fitted track, #theta_{scatt}}{measured with both Scifi and HCAL};3D scattering angle, #theta_{scatt} [degrees];Counts'
            self.hists[name] = ROOT.TH1D(name, title, 90, -90, 90)
        if not np.isnan(gradients['x']) and not np.isnan(gradients['y']):
            dotprod = slopeX*gradients['x'] + slopeY*gradients['y'] + 1
            norm = np.sqrt(slopeX**2 + slopeY**2 + 1)*np.sqrt(gradients['x']**2 + gradients['y']**2 + 1)
            scatteringangle = np.arccos(dotprod/norm)            
            self.hists[name].Fill(scatteringangle)

    def GetScifiMedians(self):
        scifi_hits = self.tw.M.eventTree.Digi_ScifiHits

        res={}
        for hit in scifi_hits:
            detID = hit.GetDetectorID()
            station = detID//1000000
            
            self.tw.M.Scifi.GetSiPMPosition(detID, self.A,self.B)
            
            if hit.isVertical(): 
                pos = 0.5*(self.A[0] + self.B[0])
                orientation=1
            else:
                pos = 0.5*(self.A[1] + self.B[1]) 
                orientation=0

            if not station in res: res[station]={'hor':[], 'vert':[]}
            if orientation==1: res[station]['vert'].append(pos)
            else: res[station]['hor'].append(pos)
        return res









    ###########
    # Method for looking for the next time an expected small SiPM fires
    # in the numu candidate events
    ###########

    def FindExpectedSmallSiPMs(self, fDetID, qdcs, smallSiPMs):

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

        if self.tw.options.signalpartitions: # If I'm looking at mu_nu candidates
            n_partition = list(self.tw.options.signalpartitions.keys()).index(f'{str(runNr).zfill(6)}')
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

    def WriteOutRecordedTimes(self):

        filename=f'/eos/home-a/aconsnd/SWAN_projects/numuInvestigation/data/numuhits'

        with open(f'{filename}.json', 'w') as f:
            json.dump(self.data, f, indent=4)
        print(f'File saved to {filename}.json')

    def SaveClusters(self):

        filename=f'/eos/home-a/aconsnd/SWAN_projects/numuInvestigation/data/clusterInfo'

        with open(f'{filename}.json', "w") as f:
            json.dump(self.clusters, f, indent=4)

        print(f'File saved to {filename}.json')                

    def SaveScifiHits(self):

        filename=f'/eos/home-a/aconsnd/SWAN_projects/numuInvestigation/data/ScifiHits'

        with open(f'{filename}.json', 'w') as f:
            json.dump(self.scifidata, f, indent=4)
        print(f'File saved to {filename}.json')                

    def GetOutFileName(self):
        if not self.simulation: 
            d = f'{self.outpath}splitfiles/run{self.runNr}/showerprofiles/'
            outfilename_noExt=f'showerprofiles_{self.tw.options.nStart}'
        else: d = f'{self.outpath}showerprofiles/'
        
        # Conditions for parsing the filenames of different simulation types
        if self.simulation and self.tw.options.simMode == 'muonDIS':
            preamble, key1, key2, end = self.tw.options.fname.split('_')
            outfilename_noExt=f'showerprofiles_{key1}_{key2}'

        elif self.simulation and self.tw.options.simMode == 'neutrino':
            simulation_details, filename = self.tw.options.fname.split('/')
            key=filename.split('_')[1]
            d+=simulation_details+'/'
            outfilename_noExt=f'showerprofiles_{key}'

        elif self.simulation and self.tw.options.simMode == 'neutralhadron':
            particle_type, Emin, Emax, key = self.tw.options.fname.split('_')[3:7]
            outfilename_noExt=f'showerprofiles_{particle_type}_{Emin}_{Emax}_{key}'

        elif self.simulation and self.tw.options.simMode == 'passingmuon':
            keys=self.tw.options.fname.split('_')[1:3]
            outfilename_noExt=f'showerprofiles_{keys[0]}_{keys[1]}'
        
        else: 
            print('fucked it bud')
            return
        return d, outfilename_noExt

    def WriteOutHistograms(self):
        
        if self.simulation: 
            self.WriteOutSimulationDetails()

        d, outfilename_noExt = self.GetOutFileName()

        if os.path.exists(d+outfilename_noExt+'.root'): outfile=ROOT.TFile.Open(d+outfilename_noExt, 'recreate')
        else: 
            path_obj=Path(d)
            path_obj.mkdir(parents=True, exist_ok=True)
            outfile=ROOT.TFile.Open(d+outfilename_noExt+'.root', 'create') 

        for hname in self.hists:

            if hname in ['polar', 'azimuthal','thetavphi', 'dXmu-polar', 'dXmu-azimuthal', 'scatteringangle']:
                outfile.WriteObject(self.hists[hname], hname)
                continue

            key = hname.split('-')[0]
            # if key=='lambda_x_maxQDCvQDC': print('HELLO THERE')
            hist=self.hists[hname]

            if not hasattr(outfile, key): 
                folder=outfile.mkdir(key)
            else: folder=outfile.Get(key)
            
            folder.cd()
            hist.Write(hname, 2) # The 2 means it will overwrite a hist of the same name  

        outfile.Close()
        print(f'{len(self.hists)} histograms saved to {d+outfilename_noExt}.root')        