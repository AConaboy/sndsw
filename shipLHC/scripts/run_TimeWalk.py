#!/usr/bin/env python
import ROOT,os,sys,subprocess,atexit,time
from XRootD import client
from XRootD.client.flags import DirListFlags, OpenFlags, MkDirFlags, QueryCode
import Monitor, SndlhcMuonReco, TimeWalk

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
parser.add_argument('-p', '--path', dest='path', help='path', type=str)
parser.add_argument("-P", "--partition", dest="partition", help="partition of data", type=int,required=False,default=-1)
parser.add_argument("-d", "--debug", dest="debug", help="debug", type=int, default=False)
parser.add_argument("-cpp", "--convRawCPP", action='store_true', dest="FairTask_convRaw", help="convert raw data using ConvRawData FairTask", default=False)
parser.add_argument( "--withCalibration", action='store_true', dest="makeCalibration", help="make QDC and TDC calibration, not taking from raw data", default=False)
parser.add_argument("--customEventChain", dest="customEventChain", action="store_true")

parser.add_argument("-f", "--inputFile", dest="fname", help="file name for MC", type=str,default=None,required=False)
parser.add_argument("-g", "--geoFile", dest="geoFile", help="geofile", required=False,default=False)
parser.add_argument("-b", "--heartBeat", dest="heartBeat", help="heart beat", default=10000,type=int)
parser.add_argument("-c", "--command", dest="command", help="command", default="")
parser.add_argument("-n", "--nEvents", dest="nEvents", help="number of events", default=-1,type=int)
parser.add_argument("-s", "--nStart", dest="nStart", help="first event", default=0,type=int)
parser.add_argument("-t", "--trackType", dest="trackType", help="DS or Scifi", default="DS")
parser.add_argument("--CorrectionType", dest="CorrectionType", help="Type of polynomial function or log function", default=5, type=int, required=False)
parser.add_argument("--Task", dest="Task", help="TimeWalk or SelectionCriteria", default="TimeWalk")
parser.add_argument("--nStations", dest="nStations", help="How many DS planes are used in the DS track fit", type=int, default=3)
parser.add_argument("--TWCorrectionRun", dest="TWCorrectionRun", help="Select what run to take TW correction parameters from. By default it is the same as the data", type=int)
parser.add_argument("--AlignmentRun", dest="AlignmentRun", help="AlignmentRun", type=int)
parser.add_argument('-D', '--datalocation', dest='datalocation', type=str, default='physics')
parser.add_argument('--state', dest='state', type=str, default='uncorrected')
parser.add_argument('--referencesystem', dest='referencesystem', type=int, default=3)

# Cuts
parser.add_argument('--OneHitPerSystem', dest='OneHitPerSystem', action='store_true')
parser.add_argument('--SlopesCut', dest='SlopesCut', action='store_true')
parser.add_argument('--nSiPMsCut', dest='nSiPMsCut', action='store_true')
parser.add_argument('--XT', dest='XT', action='store_true')

parser.add_argument('--afswork', dest='afswork', type=str, default='/afs/cern.ch/work/a/aconsnd/Timing')
parser.add_argument('--afsuser', dest='afsuser', type=str, default='/afs/cern.ch/user/a/aconsnd/Timing')
parser.add_argument('--eosH8', dest='eosH8', type=str, default='/eos/experiment/sndlhc/convertedData/commissioning/TB_H8_october/')
parser.add_argument('--eosTI18', dest='eosTI18', type=str, default='/eos/experiment/sndlhc/convertedData/commissioning/TI18/')
parser.add_argument('--mode', dest='mode', type=str, default='zeroth')
parser.add_argument('-C', '--HTCondor', dest='HTCondor', help='int (0/1), is on HTCondor?', default=0, type=int, required=False)
parser.add_argument('--numusignalevents', dest='numusignalevents', action='store_true')
parser.add_argument('--signalpartitions', dest='signalpartitions', required=False)
parser.add_argument('--numuStudy', dest='numuStudy', action='store_true')
parser.add_argument('--simulation', dest='simulation', action='store_true')
parser.add_argument('--simEngine', dest='simEngine', type=str)

parser.add_argument("--ScifiNbinsRes", dest="ScifiNbinsRes", default=100)
parser.add_argument("--Scifixmin", dest="Scifixmin", default=-2000.)
parser.add_argument("--ScifialignPar", dest="ScifialignPar", default=False)
parser.add_argument("--ScifiResUnbiased", dest="ScifiResUnbiased", default=False)
parser.add_argument("--Mufixmin", dest="Mufixmin", default=-10.)
parser.add_argument("--chi2xpred_zpos", dest="chi2xpred_zpos", help="z-position for plotting DS track red-chi2 values", type=int, default=0)
parser.add_argument("--WriteOutTrackInfo", dest="WriteOutTrackInfo", type=int, default=0)
parser.add_argument('--numbering', dest='numbering', type=str, default='systemPCB')

parser.add_argument("--goodEvents", dest="goodEvents", action='store_true',default=False)
parser.add_argument("--withTrack", dest="withTrack", action='store_true',default=False)
parser.add_argument("--nTracks", dest="nTracks",default=0,type=int)
parser.add_argument("--save", dest="save", action='store_true',default=False)
parser.add_argument("--interactive", dest="interactive", action='store_true',default=False)

parser.add_argument("--postScale", dest="postScale",help="post scale events, 1..10..100", default=-1,type=int)

options = parser.parse_args()

# if no geofile given, use defaults according to run number

if options.runNumber < 0  and not options.geoFile and not options.numusignalevents: 
    print('No run number given and no geoFile. Do not know what to do. Exit.')
    exit()
if not options.geoFile:
    if options.runNumber < 4575:
        options.geoFile =  "geofile_sndlhc_TI18_V3_08August2022.root"
    elif options.runNumber < 4855:
        options.geoFile =  "geofile_sndlhc_TI18_V5_14August2022.root"
    elif options.runNumber < 5172:
        options.geoFile =  "geofile_sndlhc_TI18_V6_08October2022.root"
    elif options.runNumber < 5485:
        options.geoFile =  "geofile_sndlhc_TI18_V7_22November2022.root"
    else:
        options.geoFile =  "geofile_sndlhc_TI18_V1_2023.root"
# to be extended for future new alignments.

# works only for runs on EOS
# if not options.numusignalevents:
if not options.server.find('eos')<0 and not options.simulation:
    if options.path.find('2023')!=-1:
        rawDataPath='/eos/experiment/sndlhc/raw_data/physics/2023/'
    elif options.path.find('2022')!=-1:
        rawDataPath='/eos/experiment/sndlhc/raw_data/physics/2022/'
    else:
        rawDataPath='/eos/experiment/sndlhc/raw_data/commissioning/TI18/data/'
    options.rawDataPath=rawDataPath
    print(f'rawDataPath: {rawDataPath}')
    runDir=rawDataPath+'run_'+str(options.runNumber).zfill(6)
    jname = "run_timestamps.json"
    dirlist  = str( subprocess.check_output("xrdfs "+os.environ['EOSSHIP']+" ls "+runDir,shell=True) ) 
    if jname in dirlist:
       with client.File() as f:          
          f.open(options.server+runDir+"/run_timestamps.json")
          status, jsonStr = f.read()
       exec("date = "+jsonStr.decode())
       options.startTime = date['start_time'].replace('Z','')
       if 'stop_time' in date:
           options.startTime += " - "+ date['stop_time'].replace('Z','')

# prepare tasks:
FairTasks = []
houghTransform = False # under construction, not yet tested
if houghTransform:
    muon_reco_task = SndlhcMuonReco.MuonReco()
    muon_reco_task.SetName("houghTransform")
    FairTasks.append(muon_reco_task)
else:
    import SndlhcTracking
    trackTask = SndlhcTracking.Tracking() 
    trackTask.SetName('simpleTracking')
    FairTasks.append(trackTask)

M = Monitor.Monitoring(options,FairTasks)
monitorTasks = {}

if options.nEvents < 0 :   options.nEvents = M.GetEntries()
if options.postScale==0 and options.nEvents>5E7: options.postScale = 100
if options.postScale==0 and options.nEvents>5E6: options.postScale = 10

if options.numusignalevents:
    import numusignals 
    numu = numusignals(options)
    numu.InvestigateSignalEvents()


if options.Task=='TimeWalk':
    if not options.mode:
        print('=='*20+f'\nNo mode given for time walk task. Give mode as zeroth, tof or tw.\n'+'=='*20)
        pyExit()
    if options.numusignalevents: options.mode='numusignalevents'
    monitorTasks['TimeWalk'] = TimeWalk.TimeWalk(options, M) 

# else:
if not options.auto:   # default online/offline mode
    # trackIDdict={'DS':3, 'Scifi':1}
    start, nEvents=options.nStart, options.nEvents
    deciles=[i/10 for i in range(11)]
    for n in range(start,start+nEvents):
        event = M.GetEvent(n)

        if ((n-start)/nEvents) in deciles:
            progstr=f'Progress: {n}/{start+nEvents}, {int(100*(n-start)/nEvents)}%'
            print('='*len(progstr)+'\n')
            print(progstr, '\n')

        for m in monitorTasks:
            monitorTasks[m].ExecuteEvent(M.eventTree)
    if 'SelectionCriteria' in monitorTasks:
        if not options.debug: 
            monitorTasks['SelectionCriteria'].WriteOutHistograms()
            if options.WriteOutTrackInfo: monitorTasks['SelectionCriteria'].SaveTrackInfos()
    if 'TimeWalk' in monitorTasks:
        if not options.debug: monitorTasks['TimeWalk'].WriteOutHistograms()
