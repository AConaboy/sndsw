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
from args_config import add_arguments

parser = ArgumentParser()

add_arguments(parser)

options = parser.parse_args()
# if no geofile given, use defaults according to run number

if options.runNumber < 0 and not options.geoFile and not options.numusignalevents: 
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
if 'TimeWalk' in monitorTasks:
    if not options.debug: monitorTasks['TimeWalk'].WriteOutHistograms()
