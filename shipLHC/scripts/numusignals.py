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

parser = ArgumentParser()

parser.add_argument("-A", "--auto", dest="auto", help="run in auto mode online monitoring",default=False,action='store_true')
parser.add_argument("--Nupdate", dest="Nupdate", help="frequence of updating online plots",default=100,type=int)
parser.add_argument("--Nlast",      dest="Nlast", help="last N events to analyze on file",default=10,type=int)
parser.add_argument("--Cosmics",      dest="cosmics", help="use default data stream if no beam",action='store_true',default=False)
parser.add_argument("--sudo", dest="sudo", help="update files on EOS",default=False,action='store_true')

parser.add_argument("-M", "--online", dest="online", help="online mode",default=False,action='store_true')
parser.add_argument("--batch", dest="batch", help="batch mode",default=False,action='store_true')
parser.add_argument("--server", dest="server", help="xrootd server",default=os.environ["EOSSHIP"])
parser.add_argument("-r", "--runNumber", dest="runNumber", help="run number", type=int,default=-1)
parser.add_argument('-p', '--path', dest='path', help='path', type=str, default='/eos/experiment/sndlhc/convertedData/physics/2022/', required=False)
parser.add_argument("-P", "--partition", dest="partition", help="partition of data", type=int,required=False,default=-1)
parser.add_argument("-d", "--debug", dest="debug", help="debug", type=int, default=False)
parser.add_argument("-cpp", "--convRawCPP", action='store_true', dest="FairTask_convRaw", help="convert raw data using ConvRawData FairTask", default=False)
parser.add_argument( "--withCalibration", action='store_true', dest="makeCalibration", help="make QDC and TDC calibration, not taking from raw data", default=False)

parser.add_argument("-f", "--inputFile", dest="fname", help="file name for MC", type=str,default=None,required=False)
parser.add_argument("-g", "--geoFile", dest="geoFile", help="geofile", required=False,default="geofile_sndlhc_TI18_V3_08August2022.root")
parser.add_argument("-b", "--heartBeat", dest="heartBeat", help="heart beat", default=10000,type=int)
parser.add_argument("-c", "--command", dest="command", help="command", default="")
parser.add_argument("-n", "--nEvents", dest="nEvents", help="number of events", default=-1,type=int)
parser.add_argument("-s", "--nStart", dest="nStart", help="first event", default=0,type=int)
parser.add_argument("-t", "--trackType", dest="trackType", help="DS or Scifi", default="DS")
parser.add_argument("--CorrectionType", dest="CorrectionType", help="Type of polynomial function or log function", default=5, type=int, required=False)
parser.add_argument("--Task", dest="Task", help="TimeWalk or SelectionCriteria", default="TimeWalk")
parser.add_argument("--numuevents", action='store_true', default=True)
parser.add_argument("--nStations", dest="nStations", help="How many DS planes are used in the DS track fit", type=int, default=3)
parser.add_argument("--TWCorrectionRun", dest="TWCorrectionRun", help="Select what run to take TW correction parameters from. By default it is the same as the data", type=int)
parser.add_argument("--AlignmentRun", dest="AlignmentRun", help="AlignmentRun", type=int)
parser.add_argument('-D', '--datalocation', dest='datalocation', type=str, default='physics')
parser.add_argument('--state', dest='state', type=str, default='uncorrected')

# Cuts
parser.add_argument('--OneHitPerSystem', dest='OneHitPerSystem', action='store_true')
parser.add_argument('--SlopesCut', dest='SlopesCut', action='store_true')
parser.add_argument('--nSiPMsCut', dest='nSiPMsCut', action='store_true')
parser.add_argument('--CrossTalk', dest='CrossTalk', action='store_true')

parser.add_argument('--afswork', dest='afswork', type=str, default='/afs/cern.ch/work/a/aconsnd/Timing')
parser.add_argument('--afsuser', dest='afsuser', type=str, default='/afs/cern.ch/work/a/aconsnd/Timing')
parser.add_argument('--eosH8', dest='eosH8', type=str, default='/eos/experiment/sndlhc/convertedData/commissioning/TB_H8_october/')
parser.add_argument('--eosTI18', dest='eosTI18', type=str, default='/eos/experiment/sndlhc/convertedData/commissioning/TI18/')
parser.add_argument('--mode', dest='mode', type=str, default='showerprofiles')
parser.add_argument('-C', '--HTCondor', dest='HTCondor', action='store_true')
parser.add_argument('--numusignalevents', dest='numusignalevents', action='store_true')
parser.add_argument('--numuStudy', dest='numuStudy', action='store_true', default=True)
# options.numuStudy=True

# Investigate small SiPMs
parser.add_argument('--SmallSiPMcheck', dest='SmallSiPMcheck', action='store_true')
# Run showerprofiles.ShowerDirection with a cut on IsTrackRelated
parser.add_argument('--dycut', dest='dycut', action='store_true')
# Run showerprofiles.ShowerDirection without considering the exact bar that the DS track extrapolates to
parser.add_argument('--notDSbar', dest='notDSbar', action='store_true')
parser.add_argument('--SiPMmediantimeCut', dest='SiPMmediantimeCut', action='store_true')
parser.add_argument('--SiPMtimeCut', dest='SiPMtimeCut', action='store_true')
# Update geofile, takes longer
parser.add_argument('--updateGeoFile', dest='updateGeoFile', action='store_true')

parser.add_argument("--ScifiNbinsRes", dest="ScifiNbinsRes", default=100)
parser.add_argument("--Scifixmin", dest="Scifixmin", default=-2000.)
parser.add_argument("--ScifialignPar", dest="ScifialignPar", default=False)
parser.add_argument("--ScifiResUnbiased", dest="ScifiResUnbiased", default=False)
parser.add_argument("--Mufixmin", dest="Mufixmin", default=-10.)
parser.add_argument("--chi2xpred_zpos", dest="chi2xpred_zpos", help="z-position for plotting DS track red-chi2 values", type=int, default=0)
parser.add_argument("--WriteOutTrackInfo", dest="WriteOutTrackInfo", type=int, default=0)
parser.add_argument('--numbering', dest='numbering', type=str, default='systemPCB')
parser.add_argument('--referencesystem', dest='referencesystem', type=int, default=3)

parser.add_argument("--goodEvents", dest="goodEvents", action='store_true',default=False)
parser.add_argument("--withTrack", dest="withTrack", action='store_true',default=False)
parser.add_argument("--nTracks", dest="nTracks",default=0,type=int)
parser.add_argument("--save", dest="save", action='store_true',default=False)
parser.add_argument("--interactive", dest="interactive", action='store_true',default=False)
parser.add_argument("--postScale", dest="postScale",help="post scale events, 1..10..100", default=-1,type=int)
parser.add_argument("--load-hists", dest="load_hists", action='store_true')

options = parser.parse_args()

class Numusignaleventtiming(object):

    def __init__(self, options=options):
       
        self.colours={0:ROOT.kRed, 1:ROOT.kBlack, 2:ROOT.kMagenta, 3: ROOT.kGreen}
        self.options=options
        self.afswork=f'{options.afswork}-physics2022/'
        self.hists={}
        # self.muAna = Analysis.Analysis(options)
        self.numucandidatescale=5
        self.legend=ROOT.TLegend(0.14, 0.60, 0.40, 0.85)
        self.planelegends={i:ROOT.TLegend(0.14, 0.60, 0.40, 0.85) for i in range(5)}
        if options.load_hists:
            self.LoadSignalHists()
            self.LoadPassingMuonHists()
            self.FillLegend()
            # self.MakeSignalComparisonCanvases()
            # self.WriteOutSignalComparisonCanvases()

        else:
            self.MakeSignalPartitions()
        options.numuStudy=True

    def MakeSignalPartitions(self):
        numusignalevent_filepath = '/afs/cern.ch/work/a/aconsnd/numusignalevents.csv'
        with open(numusignalevent_filepath, 'r') as f:
            reader=csv.reader(f)
            nu_mu_data=[r for r in reader]
        self.nu_mu_events={int(x[0]):(int(x[1]),int(x[2])) for x in nu_mu_data}
        # self.muAna.Get_numuevents()

        self.signalpartitions={}
        self.eventChain=ROOT.TChain("rawConv")
        
        for runNr in self.nu_mu_events:
            eventNumber=self.nu_mu_events[runNr][0]

            runNumber = str(runNr).zfill(6)
            partition = int(eventNumber // 1E6)
            dirlist  = os.listdir(options.path+"run_"+runNumber)
            partitionfile=f'{options.path}run_{runNumber}/sndsw_raw-{str(partition).zfill(4)}.root'        
            self.signalpartitions[runNumber]=str(partition).zfill(4)
            self.eventChain.Add(partitionfile)

        options.signalpartitions = self.signalpartitions
        options.customEventChain=self.eventChain

        # Set initial geofile to instance Monitor
        options.runNumber = list(self.nu_mu_events.keys())[0]

        FairTasks=[]
        trackTask = SndlhcTracking.Tracking()
        trackTask.SetName('simpleTracking')
        FairTasks.append(trackTask)  

        self.M = Monitor.Monitoring(options, FairTasks)

        self.monitorTasks = {}

        if options.Task=='TimeWalk':
            options.mode='showerprofiles'
            self.monitorTasks['TimeWalk'] = TimeWalk.TimeWalk(options, self.M) 
              
        self.tw = self.monitorTasks['TimeWalk']
        self.muAna = self.tw.muAna
        
        self.hists=self.tw.hists

    def InvestigateSignalEvents(self):

        runs=self.nu_mu_events.keys()

        for runNr in runs:
            self.InvestigateEvent(runNr)
        self.tw.sp.WriteOutRecordedTimes()
        self.MakeAngularPlots()
        
    def InvestigateEvent(self, runNr):
        evt_number=self.GetSignalEventNumber(runNr)
        self.M.GetEvent(evt_number) # Runs tracking task
        
        eventheader = self.M.eventTree.EventHeader
        runId = eventheader.GetRunId()
        
        # print(f'RunNr: {runNr}, runId: {runId}')
        # print(f'signal event in M.eventTree: {evt_number}, M.EventNumber: {self.M.EventNumber}, event number in partition: {eventheader.GetEventNumber()}')

        if options.updateGeoFile:
            if runNr < 4575:     options.geoFile =  "geofile_sndlhc_TI18_V3_08August2022.root"
            elif runNr < 4855:   options.geoFile =  "geofile_sndlhc_TI18_V5_14August2022.root"
            elif runNr < 5172:   options.geoFile =  "geofile_sndlhc_TI18_V6_08October2022.root"
            elif runNr < 5485:   options.geoFile =  "geofile_sndlhc_TI18_V7_22November2022.root"
            else:                options.geoFile =  "geofile_sndlhc_TI18_V1_2023.root"            

            self.M.snd_geo = SndlhcGeo.GeoInterface(options.path+options.geoFile)
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
        print('Alignment:', alignment)
        self.muAna.MakeAlignmentParameterDict(alignment) # Load appropriate alignment parameters into muAna
        self.muAna.MakeTWCorrectionDict(alignment)

        for m in self.monitorTasks:
            self.monitorTasks[m].ExecuteEvent(self.M.eventTree)
            
    def GetSignalEventNumber(self, runNr):
        n_partition = list(self.nu_mu_events.keys()).index(runNr)
        evt_number = int(n_partition * 1e6 + self.nu_mu_events[runNr][0] % 1e6)            
        return evt_number
    
    def testing(self):
        self.fired_detIDs={}
        for run in self.nu_mu_events:
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
        'dx':[], 'ybarycentre':[], 'dy':[], 'xEx':[], 'yEx':[], 'interaction wall':[]}
        
        """
        sp.barycentres = {
                            runNr:
                                {
                                    plane: [[x_barycentre, y_barycentre], [xEx, yEx]]
                                }
                        }
        """

        for runNr in self.tw.sp.barycentres.keys():
            barycentres = self.tw.sp.barycentres[runNr] # { plane: [ [(x_barycentre,dxbc), (y_barycentre,dxbc)], [xEx, yEx]] }
            
            planes = [i for i in barycentres.keys() if len(barycentres[i])!=0]
            
            xbcs = [barycentres[i][0][0] for i in planes]
            xbarycentres, dxs = zip(*xbcs)
            
            ybcs = [barycentres[i][0][1] for i in planes]
            ybarycentres, dys = zip(*ybcs)

            extrapolations = [barycentres[i][1] for i in planes]
            xExs, yExs = zip(*extrapolations)

            interactionwall = self.nu_mu_events[runNr][1]

            self.data['runNr'].append(runNr)
            self.data['planes'].append(planes)
            self.data['xbarycentre'].append(xbarycentres)
            self.data['dx'].append(dxs) # uncertainty on xbarycentre
            self.data['ybarycentre'].append(ybarycentres)
            self.data['dy'].append(dys) # uncertainty on xbarycentre
            self.data['xEx'].append(xExs) # DS track extrapolated to plane 
            self.data['yEx'].append(yExs) # DS track extrapolated to plane
            self.data['interaction wall'].append(interactionwall)

        self.df=pd.DataFrame(self.data)

        filename='/eos/home-a/aconsnd/SWAN_projects/numuInvestigation/data/barycentres'
        if self.options.notDSbar: filename+='-notDSbar'
        elif self.options.dycut: filename+='-dycut'
        
        if self.options.SiPMmediantimeCut: filename+='-SiPMmediantimeCut'
        if self.options.SiPMtimeCut: filename+='-SiPMtimeCut'
        self.df.to_csv(f'{filename}.csv')
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