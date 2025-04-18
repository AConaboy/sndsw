#!/usr/bin/env python
import ROOT,os,sys, csv, subprocess,atexit,time
from XRootD import client
from XRootD.client.flags import DirListFlags, OpenFlags, MkDirFlags, QueryCode
import Monitor, SndlhcMuonReco, SndlhcTracking, SndlhcGeo, TimeWalk
import AnalysisFunctions as Analysis 
import matplotlib.pyplot as plt 
import matplotlib.gridspec as gridspec

ROOT.gROOT.SetBatch(True)

def pyExit():
    print("Make suicide until solution found for freezing")
    os.system('kill '+str(os.getpid()))
atexit.register(pyExit)

from argparse import ArgumentParser
from args_config import add_arguments

parser = ArgumentParser()
parser.add_argument('--analysisName', dest='analysisName', type=str, default='numu')
parser.add_argument('--load_hists', dest='load_hists', action='store_true')
parser.add_argument('--updateGeoFile', dest='updateGeoFile', action='store_true')
add_arguments(parser)

options = parser.parse_args()

if options.analysisName=='numu': options.numuStudy=True 
elif options.analysisName=='nue': options.nueStudy=True 

class Numusignaleventtiming(object):

    def __init__(self, options=options):
       
        self.colours={0:ROOT.kRed, 1:ROOT.kBlack, 2:ROOT.kMagenta, 3: ROOT.kGreen}
        self.options=options
        self.afswork=f'{options.afswork}-physics2022/'
        self.hists={}
        # self.muAna = Analysis.Analysis(options)
        self.legend=ROOT.TLegend(0.14, 0.60, 0.40, 0.85)
        self.planelegends={i:ROOT.TLegend(0.14, 0.60, 0.40, 0.85) for i in range(5)}

        if options.load_hists:
            self.LoadSignalHists()
            self.LoadPassingMuonHists()
            self.FillLegend()

        else:
            self.MakeSignalPartitions()

    def MakeSignalPartitions(self):
        
        self.signal_events={}

        if options.analysisName=='numu': 
            signalevent_filepath = '/afs/cern.ch/work/a/aconsnd/numusignalevents.csv'
        
            with open(signalevent_filepath, 'r') as f:
                reader=csv.reader(f)
                for idx,x in enumerate(reader):
                    if idx==0: continue
                    self.signal_events[int(x[0])] = [int(x[1]), int(x[2])] + [float(i) for i in x[3:]]

        elif options.analysisName=='nue':
            signalevent_filepath = '/afs/cern.ch/work/a/aconsnd/nuesignalevents.csv'

            with open(signalevent_filepath, 'r') as f:
                reader=csv.reader(f)
                for idx,x in enumerate(reader):
                    if idx==0: continue
                    self.signal_events[int(x[0])] = [int(x[1])]     

        self.signalpartitions={}
        self.eventChain=ROOT.TChain("rawConv")
        # self.eventChain.SetBranchStatus("EventHeader", 0) # To be able to run on old MC.
        # self.eventChain.SetBranchStatus("EventHeader.", 0) # To be able to run on old MC. Bizarre....?
        
        for runNr in self.signal_events:
            year = self.GetRunYear(runNr)
            path = f'/eos/experiment/sndlhc/convertedData/physics/{year}'
            eventNumber=self.signal_events[runNr][0]
            runNumber = str(runNr).zfill(6)
            partition = int(eventNumber // 1E6)

            partitionfile=f'{path}/run_{runNumber}/sndsw_raw-{str(partition).zfill(4)}.root'  
            print(partitionfile)
            if not os.path.exists(partitionfile):
                print(f'No file: {partitionfile}')
                continue
            self.signalpartitions[runNumber]=str(partition).zfill(4)
            self.eventChain.Add(partitionfile)

        options.signalpartitions = self.signalpartitions
        options.customEventChain = self.eventChain

        # Set initial geofile to instance Monitor
        options.runNumber = list(self.signal_events.keys())[0]
        options.geoFile = 'geofile_sndlhc_TI18.root'
        start_year = self.GetRunYear(options.runNumber)
        options.path = f'/eos/experiment/sndlhc/convertedData/physics/{start_year}/'

        FairTasks=[]
        trackTask = SndlhcTracking.Tracking()
        trackTask.SetName('simpleTracking')
        FairTasks.append(trackTask)

        self.M = Monitor.Monitoring(options, FairTasks)

        self.monitorTasks = {}

        if options.Task=='TimeWalk':

            if options.analysisName=='numu':options.mode='showerprofiles'
            if options.analysisName=='nue':options.mode='extendedreconstruction'
            self.monitorTasks['TimeWalk'] = TimeWalk.TimeWalk(options, self.M) 
            
        self.tw = self.monitorTasks['TimeWalk']
        self.muAna = self.tw.muAna
        
        self.hists=self.tw.hists

    def GetRunYear(self, runNr):
        if isinstance(runNr, str): runNr=int(runNr)

        if runNr < 5485: year='2022'
        elif 5485 <= runNr < 7656 : year='2023'
        else: year='2024'

        return year

    def InvestigateSignalEvents(self):

        runs=self.signal_events.keys()
        self.muAna.MakeTWCorrectionDict()

        for runNr in runs:
            if runNr==4752: continue
            tmp = self.InvestigateEvent(runNr)
            if not all(self.tw.sp.pass_cuts.values()):
                # print(f'Not all cuts are passing! Must be wrong event for run {runNr}!')
                print(self.tw.sp.pass_cuts)
            if tmp==-999: return
        self.tw.sp.WriteOutRecordedTimes()
        self.tw.sp.SaveScifiHits()
        # self.MakeAngularPlots()
        
    def InvestigateEvent(self, runNr):
        evt_number=self.GetSignalEventNumber(runNr)
        self.M.GetEvent(evt_number) # Runs tracking task

        eventheader = self.M.eventTree.EventHeader
        eventheader_eventnumber = self.M.eventTree.EventHeader.GetEventNumber()
        runId = eventheader.GetRunId()
        if not runId==runNr:
            print(f'FUCKED. runId={runId}, runNr={runNr}')
            return

        if eventheader_eventnumber != self.signal_events[runNr][0]:
            print(f'Event number in run not equal to event number in event header!')        
            print(f'run event number: {self.signal_events[runNr][0]}, event header number: {eventheader_eventnumber}')
            print(f'==='*5)

        if options.updateGeoFile:
            if runNr < 4575:     options.geoFile =  "geofile_sndlhc_TI18_V3_08August2022.root"
            elif runNr < 4855:   options.geoFile =  "geofile_sndlhc_TI18_V5_14August2022.root"
            elif runNr < 5172:   options.geoFile =  "geofile_sndlhc_TI18_V6_08October2022.root"
            elif runNr < 5485:   options.geoFile =  "geofile_sndlhc_TI18_V7_22November2022.root"
            else:                options.geoFile =  "geofile_sndlhc_TI18_V1_2023.root"            

            year = self.GetRunYear(runNr)
            path = f'/eos/experiment/sndlhc/convertedData/physics/{year}'
            self.M.snd_geo = SndlhcGeo.GeoInterface(path+options.geoFile)
            self.M.MuFilter = self.M.snd_geo.modules['MuFilter']
            self.M.Scifi       = self.M.snd_geo.modules['Scifi']
            self.M.zPos = self.M.getAverageZpositions()
        
        # print(f'Investigating run {runNr}, run event {evt_number}')
        firedUSDetIds = [i.GetDetectorID() for i in self.M.eventTree.Digi_MuFilterHits if i.GetDetectorID()//10000==2]
        nVeto = len([i.GetDetectorID() for i in self.M.eventTree.Digi_MuFilterHits if i.GetDetectorID()//10000==1])
        
        # Sanity check that correct event is loaded
        if nVeto != 0: 
            print(f'Veto in this event. Cannot be correct event')

        # Get time alignment type of signal event
        alignment = self.muAna.GetTimeAlignmentType(runId)
        print(f'run: {runId}, alignment: {alignment}')
        self.tw.timealignment=alignment
        self.muAna.timealignment=alignment
        self.muAna.MakeAlignmentParameterDict() # Load appropriate alignment parameters into muAna
        self.muAna.MakeTimingCovarianceDict()

        for m in self.monitorTasks:
            self.monitorTasks[m].ExecuteEvent(self.M.eventTree)
            
    def GetSignalEventNumber(self, runNr):
        # This gets the event number in the custom chain! 
        n_partition = list(self.signal_events.keys()).index(runNr)
        evt_number = int(n_partition * 1e6 + self.signal_events[runNr][0] % 1e6)
        return evt_number
    
    def testing(self):
        self.fired_detIDs={}
        for run in self.signal_events:
            signalevent = self.GetSignalEventNumber(run)
            self.M.GetEvent(signalevent)

            hits = self.M.eventTree.Digi_MuFilterHits
            self.fired_detIDs[run] = [i.GetDetectorID() for i in hits if i.GetDetectorID()//10000!=3]
            
    def LoadSignalHists(self):
        fname = f'{self.afswork}Results/SignalComparisonPlots.root'
        f=ROOT.TFile.Open(fname, 'read')

        for idx,dirname in enumerate(('averagetime', 'sidetime', 'deltatime')):

            tdir = f.Get(dirname)
            tdir.cd()

            for jdx, key in enumerate(tdir.GetListOfKeys()):
                hist=tdir.Get(key.GetName()).Clone()
                hist.SetDirectory(ROOT.gROOT)
                hist.SetLineColor(ROOT.kRed)
                hist.SetName(f'signal-{hist.GetName()}')
                self.hists[hist.GetName()]=hist

        self.planesignalhists={}
        self.planesignalstacks={}
        for histname in self.hists:
            tmp, detID, state=histname.split('_')
            key=tmp.split('-')[1]
            if key=='sidetime': 
                detID, side = detID.split('-')
                key = f'sidetime-{side}'
            s,p,b = self.muAna.parseDetID(detID)

            if not p in self.planesignalhists: 
                self.planesignalhists[p]={}
                self.planesignalstacks[p]={}
            if not key in self.planesignalhists[p]: 
                self.planesignalhists[p][key]=[]
            self.planesignalhists[p][key].append(self.hists[histname])

        for p in self.planesignalhists:
            for key in self.planesignalhists[p]:
                for idx, hist in enumerate(self.planesignalhists[p][key]): hist.SetLineColor(self.colours[idx]) 

    def CheckEventNumber(self, runNr):
        evt_number = self.GetSignalEventNumber(runNr)
        self.M.eventTree.GetEvent(evt_number)

        eventheader_eventnumber = self.M.eventTree.EventHeader.GetEventNumber()
        event_inRun = self.signal_events[runNr][0]
        
        sameEvent = eventheader_eventnumber==event_inRun
        print(f'{event_inRun}, {eventheader_eventnumber}, {sameEvent}')

    def LoadPassingMuonHists(self):

        self.fired_detIDs = []

        signalhists=self.hists
        for i in signalhists:

            if len(i.split('_'))!=3: continue
            if i.lower().find('smallsipm'): continue
            key, detID, state = i.split('_')
            if len(detID.split('-')) != 1: continue
            
            if detID in self.fired_detIDs: continue
            self.fired_detIDs.append(detID)

        self.muonhists={}

        for idx, detID in enumerate(self.fired_detIDs):
            
            fname = f'{self.afswork}/rootfiles/run005408/timewalk_{detID}_0.root'
            f=ROOT.TFile.Open(fname, 'read')

            averagetime_name = f'averagetime_{detID}_aligned'
            avgtime_dir = f.Get('averagetime')
            avgtime_dir.cd()
            avgtime_hist=avgtime_dir.Get(averagetime_name).Clone()
            avgtime_hist_dt = avgtime_hist.ProjectionY(f'muons-{avgtime_hist.GetName()}-dt')
            avgtime_hist_dt.SetDirectory(ROOT.gROOT)
            self.muonhists[avgtime_hist_dt.GetName()] = avgtime_hist_dt

            deltatime_name = f'deltatime_{detID}_aligned'
            deltatime_dir = f.Get('deltatime')
            deltatime_dir.cd()
            deltatime_hist=deltatime_dir.Get(deltatime_name).Clone()
            deltatime_hist_dt = deltatime_hist.ProjectionY(f'muons-{deltatime_hist.GetName()}-dt')
            deltatime_hist_dt.SetDirectory(ROOT.gROOT)
            self.muonhists[deltatime_hist_dt.GetName()] = avgtime_hist_dt

            avgsidetimedir = f.Get('sidetime')
            avgsidetimedir.cd()
            for side in ('left', 'right'):
                averagesidetime_name = f'sidetime_{detID}-{side}_aligned'
                averagesidetime_hist = avgsidetimedir.Get(averagesidetime_name).Clone()
                averagesidetime_hist_dt = averagesidetime_hist.ProjectionY(f'muons-{averagesidetime_hist.GetName()}-dt')
                averagesidetime_hist_dt.SetDirectory(ROOT.gROOT)
                self.muonhists[averagesidetime_hist_dt.GetName()] = averagesidetime_hist_dt

            f.Close()

        self.planemuonhists={}
        for histname in self.muonhists:
            tmp, detID, state=histname.split('_')
            key=tmp.split('-')[1]
            if key=='sidetime': 
                detID, side = detID.split('-')
                key = f'sidetime-{side}'
            s,p,b = self.muAna.parseDetID(detID)

            if not p in self.planemuonhists: self.planemuonhists[p]={}
            if not key in self.planemuonhists[p]: self.planemuonhists[p][key]=[]
            self.planemuonhists[p][key].append(self.muonhists[histname])

        for p in self.planemuonhists:
            for key in self.planemuonhists[p]:
                
                if key.find('sidetime')==0: side=key.split('-')[1]

                hists=self.planemuonhists[p][key]
                mergedhist=hists[0].Clone()
                [mergedhist.Add(i) for i in hists]
                if key.find('sidetime')==0: mergedhist.SetName(f'sidetime_plane{p}-{side}_aligned')
                else: mergedhist.SetName(f'{key}_plane{p}_aligned')
                self.planemuonhists[p][key]=mergedhist

    def FillLegend(self):
        signalhist=self.hists[list(self.hists.keys())[0]]
        muonhist=self.muonhists[list(self.muonhists.keys())[0]]

        self.legend.AddEntry(signalhist, '#nu_{#mu} data')
        self.legend.AddEntry(muonhist, '#mu data')

        for p in self.planesignalhists:
            for hist in self.planesignalhists[p]['averagetime']:
                detID = hist.GetName().split('_')[1]
                s,p,b=self.muAna.parseDetID(detID)
                self.planelegends[p].AddEntry(hist, '#nu_{#mu} data, bar '+str(b+1))

    def MakeSignalComparisonCanvases(self):
        
        self.secondaxes={}
        self.comparisoncanvases={}
        self.planecomparisoncanvases={}
        
        for idx, detID in enumerate(self.fired_detIDs):
            s,p,b = self.muAna.parseDetID(detID)

            canvasname=f'signalvmuon_comparison_{detID}'
            canvas=ROOT.TCanvas(canvasname, canvasname, 1000, 500)
            canvas.Divide(2, 2)
            self.comparisoncanvases[canvasname]=canvas

            if not p in self.planecomparisoncanvases:
                canvasname=f'signalvmuon_comparison_plane{p}'
                self.planecomparisoncanvases[canvasname]=ROOT.TCanvas(canvasname, canvasname, 1000, 500)
                self.planecomparisoncanvases[canvasname].Divide(2,2)
                for jdx, histkey in enumerate(('averagetime', 'deltatime', 'left', 'right')):
                    self.MakeCompositeHist(detID, histkey, plane=True)
            
            for jdx, histkey in enumerate(('averagetime', 'deltatime', 'left', 'right')):
                self.MakeCompositeHist(detID, histkey)

        print(f'Comparison canvases made')

    def MakeCompositeHist(self, detID, histkey, plane=None):

        padkeys={'averagetime':1, 'deltatime':2, 'left':3, 'right':4}
        print(f'Making composite hist for {detID}, {histkey}, plane:{plane}')
        if not plane: readabledetID = self.muAna.MakeHumanReadableDetID(detID)
        
        else:
            s,p,b=self.muAna.parseDetID(detID) 
            readabledetID = 'upstream, plane '+str(p+1)
        
        titlepart1='Comparison of signal times against passing muon times'
        if histkey=='averagetime':  
            if not plane: 
                signalhistname = f'signal-averagetime_{detID}_aligned'
                muonhistname = f'muons-averagetime_{detID}_aligned-dt'
                title='#splitline{'+titlepart1+'}{'+readabledetID+' average time};Bar average time [ns];Counts'
            else:
                title='#splitline{'+titlepart1+'}{'+readabledetID+' average time};Plane average time [ns];Counts'
        
        elif histkey=='deltatime':  
            if not plane: 
                signalhistname = f'signal-averagetime_{detID}_aligned'
                muonhistname = f'muons-averagetime_{detID}_aligned-dt'
                title='#splitline{'+titlepart1+'}{'+readabledetID+' average time};Plane #Delta side time [ns];Counts'
            else:
                title='#splitline{'+titlepart1+'}{'+readabledetID+' average time};Plane #Delta side time [ns];Counts'            
        
        elif histkey=='left':  
            if not plane: 
                signalhistname = f'signal-averagetime_{detID}_aligned'
                muonhistname = f'muons-averagetime_{detID}_aligned-dt'
                title='#splitline{'+titlepart1+'}{'+readabledetID+' average time};Plane left side time [ns];Counts'
            else:
                title='#splitline{'+titlepart1+'}{'+readabledetID+' average time};Plane left side time [ns];Counts'
        
        elif histkey=='right':  
            if not plane: 
                signalhistname = f'signal-averagetime_{detID}_aligned'
                muonhistname = f'muons-averagetime_{detID}_aligned-dt'
                title='#splitline{'+titlepart1+'}{'+readabledetID+' average time};Plane right side time [ns];Counts'
            else:
                title='#splitline{'+titlepart1+'}{'+readabledetID+' average time};Plane right side time [ns];Counts'

        if not plane:
            signalhist = self.hists[signalhistname]
            muonhist = self.muonhists[muonhistname]
        else:
            if histkey in ('left', 'right'):
                signalhists = self.planesignalhists[p][f'sidetime-{histkey}']
                muonhist = self.planemuonhists[p][f'sidetime-{histkey}']
            else: 
                signalhists = self.planesignalhists[p][f'{histkey}']
                muonhist = self.planemuonhists[p][f'{histkey}']                

        if not plane: pad=self.comparisoncanvases[f'signalvmuon_comparison_{detID}'].cd(padkeys[histkey])
        else: pad=self.planecomparisoncanvases[f'signalvmuon_comparison_plane{p}'].cd(padkeys[histkey])
        muonhist.SetTitle(title)
        muonhist.SetTitleSize(0.06, 'xyzt')
        muonhist.SetStats(0)
        muonhist.Draw()

        pad.Update()

        rightmax=1.1*muonhist.GetMaximum()
        scale=ROOT.gPad.GetUymax()/rightmax
        
        if plane: 
            for idx,h in enumerate(signalhists):
                h.Scale(scale/rightmax)
                h.Draw('same')

        else: signalhist.Scale(scale/rightmax)
        
        if not plane: self.legend.Draw()
        else: self.planelegends[p].Draw()

    def MakeAngularPlots(self):
        import pandas as pd

        zPositions = [self.M.zPos['MuFilter'][20+i] for i in range(5)]

        self.data={'runNr':[], 'planes':[], 'xbarycentre':[],
        'dx':[], 'dy':[], 'xEx':[], 'yEx':[], 'interaction wall':[]}
        
        """
        sp.barycentres = {
                            runNr:
                                {
                                    plane: [[x_barycentre, y_barycentre], [xEx, yEx]]
                                }
                        }
        """

        for runNr in self.tw.sp.all_barycentres.keys():
            barycentres = self.tw.sp.all_barycentres[runNr] # { plane: [ [(x_barycentre,dxbc), (y_barycentre,dxbc)], [xEx, yEx]] }
            
            planes = [i for i in barycentres.keys() if len(barycentres[i])!=0]
            
            xbcs = [barycentres[i][0][0] for i in planes]
            xbarycentres, dxs = zip(*xbcs)
            
            ybcs = [barycentres[i][0][1] for i in planes]
            ybarycentres, dys = zip(*ybcs)

            extrapolations = [barycentres[i][1] for i in planes]
            xExs, yExs = zip(*extrapolations)

            interactionwall = self.signal_events[runNr][1]

            self.data['runNr'].append(runNr)
            self.data['planes'].append(planes)
            self.data['xbarycentre'].append(xbarycentres)
            self.data['dx'].append(dxs) # uncertainty on xbarycentre
            self.data['ybarycentre'].append(ybarycentres)
            self.data['dy'].append(dys) # uncertainty on xbarycentre
            self.data['xEx'].append(xExs) # DS track extrapolated to plane 
            self.data['yEx'].append(yExs) # DS track extrapolated to plane
            self.data['interaction wall'].append(interactionwall)

        self.barycentres_df=pd.DataFrame(self.data)

        filename='/eos/home-a/aconsnd/SWAN_projects/numuInvestigation/data/barycentres'

        self.barycentres_df.to_csv(f'{filename}.csv')
        print(f'Data written to {filename}.csv')

    def WriteOutSignalHists(self):

        outfilename = f'{self.tw.outpath}Results/SignalComparisonPlots.root'
        outfile=ROOT.TFile.Open(outfilename, 'recreate')

        for histname in self.tw.hists:
            if histname in ('ExtraHitsMultiplicity', 'DISradius_US', 'US-SiPMQDC', 'frac_SiPMs', 'fractionMissingSmallFound', 'fractionMissingSmallFound-highQDC'):
                hist=self.tw.hists[histname]
                outfile.WriteObject(hist, hist.GetName(), 'kOverwrite')

            if len(histname.split('_')) == 3:

                hist=self.tw.hists[histname]
                key, detID, state=histname.split('_')
                
                if not hasattr(outfile, key): outfile.mkdir(key)
                tdir = outfile.Get(key)
                tdir.cd()

                hist.Write(hist.GetName(), 2)

            elif len(histname.split('_')) == 2:

                hist=self.tw.hists[histname]
                key, whatever = histname.split('_')
                
                if not hasattr(outfile, key): outfile.mkdir(key)
                tdir = outfile.Get(key)
                tdir.cd()
                
                hist.Write(hist.GetName(), 2)

        outfile.Close()
        print(f'Signal hists written to {outfilename}')

    def WriteOutSignalComparisonCanvases(self):

        outfilename = f'{self.afswork}Results/SignalComparisonPlots.root'
        outfile=ROOT.TFile.Open(outfilename, 'update')
        if not hasattr(outfile, 'comparison'): outfile.mkdir('comparison')
        comparisondir = outfile.Get('comparison')
        comparisondir.cd()

        for canvasname in self.comparisoncanvases:
            canvas=self.comparisoncanvases[canvasname]
            canvas.Write(canvas.GetName(), 2)

        if not hasattr(outfile, 'plane-comparison'): outfile.mkdir('plane-comparison')
        comparisondir = outfile.Get('plane-comparison')
        comparisondir.cd()

        for canvasname in self.planecomparisoncanvases:
            canvas=self.planecomparisoncanvases[canvasname]
            canvas.Write(canvas.GetName(), 2)

        outfile.Close()
        print(f'Signal comparison canvases written to {outfilename}')

if options.HTCondor==1:
    numu=Numusignaleventtiming()