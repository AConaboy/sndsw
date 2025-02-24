#!/usr/bin/env python
firstEvent = 0

import resource
def mem_monitor():
 # Getting virtual memory size 
		pid = os.getpid()
		with open(os.path.join("/proc", str(pid), "status")) as f:
				lines = f.readlines()
		_vmsize = [l for l in lines if l.startswith("VmSize")][0]
		vmsize = int(_vmsize.split()[1])
		#Getting physical memory size  
		pmsize = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
		print("memory: virtuell = %5.2F MB  physical = %5.2F MB"%(vmsize/1.0E3,pmsize/1.0E3))

import ROOT,os,sys
import shipRoot_conf
import shipunit as u
from pathlib import Path

shipRoot_conf.configure()

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("-f", "--inputFile", dest="inputFile", help="single input file", required=True)
parser.add_argument("-g", "--geoFile", dest="geoFile", help="geofile", required=True)
parser.add_argument("-p", "--path",dest="path",  help="Output directory", default='/eos/experiment/sndlhc/users/aconsnd/simulation/', required=False)
parser.add_argument("-n", "--nEvents", dest="nEvents",  type=int, help="number of events to process", default=100000)
parser.add_argument("-ts", "--thresholdScifi", dest="ts", type=float, help="threshold energy for Scifi [p.e.]", default=3.5)
parser.add_argument("-ss", "--saturationScifi", dest="ss", type=float, help="saturation energy for Scifi [p.e.]", default=104.)
parser.add_argument("-tML", "--thresholdMufiL", dest="tml", type=float, help="threshold energy for Mufi large [p.e.]", default=0.0)
parser.add_argument("-tMS", "--thresholdMufiS", dest="tms", type=float, help="threshold energy for Mufi small [p.e.]", default=0.0)
parser.add_argument("-no-cls", "--noClusterScifi", action='store_true', help="do not make Scifi clusters")
parser.add_argument("-cpp", "--digiCPP", action='store_true', dest="FairTask_digi", help="perform digitization using DigiTaskSND")
parser.add_argument("-d", "--Debug", dest="debug", action="store_true")
parser.add_argument('--simEngine', dest='simEngine', type=str)

options = parser.parse_args()
# rephrase the no-cluster flag
makeClusterScifi = not options.noClusterScifi
# -----Timer-------------
timer = ROOT.TStopwatch()
timer.Start()

# If working with Daniele's data, output path must be changed
if options.inputFile.find('dancc/MuonDIS') != -1:
	f = options.inputFile.split('/')[-1]
	file_number = options.inputFile.split('/')[-2]
	output_file = options.path + 'muonDIS/data/' + f.replace('.root', f'_{file_number}.root')
		

# Working with neutral hadrons
elif options.inputFile.find('NeutralHadrons') != -1:
	f = options.inputFile.split('/')[-1]
	
	particle_type, Emin, Emax = options.inputFile.split('/')[-4].split('_')
	
	# sanity check
	if not particle_type in ('K', 'neu'):
		print(f'Input file name not compatible!\nExiting...')
		sys.exit(1)

	file_number=options.inputFile.split('/')[-2]
	output_file = options.path+'neutralhadron/data/'+f.replace('.root', f'_{particle_type}_{Emin}_{Emax}_{file_number}.root')

elif options.inputFile.find('MuonBackground') != -1:

	f=options.inputFile.split('/')[-1]
	sim_folder=options.inputFile.split('/')[-3]
	key=options.inputFile.split('/')[-2]

	output_file = options.path+'passingmuon/data/'+f.replace('.root', f'_{sim_folder}_{key}.root')

# Working with neutrino data 
elif options.inputFile.find('/Neutrinos/Genie') != -1:
	f = options.inputFile.split('/')[-1]
	keys, file_number = options.inputFile.split('/')[-3: -1]

	output_file = options.path+'neutrino/data/'+keys+'/'+f.replace('.root', f'_{file_number}.root')
	
elif options.inputFile.find('/')==-1:
	f = options.inputFile.split('/')[-1]
	output_file = options.inputFile

output_file = output_file.replace('.root','_digi.root')  
print(f'output file: {output_file}')
# sys.exit(1)

# All data is stored on /eos
if options.inputFile.find('/eos')==0:
	if options.FairTask_digi: # flag for C++ digitisation, should be used
		options.inputFile = os.environ['EOSSHIP']+options.inputFile
	else:   
		os.system('xrdcp '+os.environ['EOSSHIP']+options.inputFile+' '+output_file)
else:
	if not options.FairTask_digi:
		os.system('cp '+options.inputFile+' '+outFile)    

# -----Create geometry----------------------------------------------
import shipLHC_conf as sndDet_conf

if options.geoFile.find('/eos')==0:
	options.geoFile = os.environ['EOSSHIP']+options.geoFile
import SndlhcGeo
snd_geo = SndlhcGeo.GeoInterface(options.geoFile)

# set digitization parameters for MuFilter
lsOfGlobals = ROOT.gROOT.GetListOfGlobals()
scifiDet = lsOfGlobals.FindObject('Scifi')
mufiDet = lsOfGlobals.FindObject('MuFilter')
mufiDet.SetConfPar("MuFilter/DsAttenuationLength",350 * u.cm)		#  values between 300 cm and 400cm observed for H6 testbeam
mufiDet.SetConfPar("MuFilter/DsTAttenuationLength",700 * u.cm)		# top readout with mirror on bottom
mufiDet.SetConfPar("MuFilter/VandUpAttenuationLength",999 * u.cm)	# no significante attenuation observed for H6 testbeam
mufiDet.SetConfPar("MuFilter/DsSiPMcalibrationS",25.*1000.)			# in MC: 1.65 keV are about 41.2 qdc
mufiDet.SetConfPar("MuFilter/VandUpSiPMcalibration",25.*1000.)
mufiDet.SetConfPar("MuFilter/VandUpSiPMcalibrationS",25.*1000.)
# mufiDet.SetConfPar("MuFilter/VandUpPropSpeed",12.5*u.cm/u.nanosecond);
mufiDet.SetConfPar("MuFilter/DsPropSpeed",14.3*u.cm/u.nanosecond)
scifiDet.SetConfPar("Scifi/nphe_min",options.ts)   # threshold
scifiDet.SetConfPar("Scifi/nphe_max",options.ss) # saturation
scifiDet.SetConfPar("Scifi/timeResol",150.*u.picosecond) # time resolution in ps
# scifiDet.SetConfPar("MuFilter/timeResol",150.*u.picosecond) # time resolution in ps, first guess

# Fair digitization task
if options.FairTask_digi:
	output_file = output_file.replace('.root','CPP.root')
	run = ROOT.FairRunAna()
	ioman = ROOT.FairRootManager.Instance()
	ioman.RegisterInputObject('Scifi', snd_geo.modules['Scifi'])
	ioman.RegisterInputObject('MuFilter', snd_geo.modules['MuFilter'])
	# Don't use FairRoot's default event header settings
	run.SetEventHeaderPersistence(False)
	
	# Set input
	fileSource = ROOT.FairFileSource(options.inputFile)
	run.SetSource(fileSource)
	# Set output
	outfile = ROOT.FairRootFileSink(output_file)
	run.SetSink(outfile)

	# Set number of events to process
	inRootFile = ROOT.TFile.Open(options.inputFile)
	inTree = inRootFile.Get('cbmsim')
	nEventsInFile = inTree.GetEntries()
	nEvents = min(nEventsInFile, options.nEvents)

	rtdb = run.GetRuntimeDb()
	DigiTask = ROOT.DigiTaskSND()
	DigiTask.withScifiClusters(makeClusterScifi)
	run.AddTask(DigiTask)
	run.Init()
	run.Run(firstEvent, nEvents)

# Digitization using python code SndlhcDigi
else:
 # import digi task
	import SndlhcDigi
	Sndlhc = SndlhcDigi.SndlhcDigi(output_file,makeClusterScifi)

	nEvents = min(Sndlhc.sTree.GetEntries(),options.nEvents)
# main loop
	for iEvent in range(firstEvent, nEvents):
		if iEvent % 50000 == 0 or options.debug:
				print('event ', iEvent, nEvents - firstEvent)
		Sndlhc.iEvent = iEvent
		rc = Sndlhc.sTree.GetEvent(iEvent)
		Sndlhc.digitize()
		if makeClusterScifi:
			 Sndlhc.clusterScifi()
 # memory monitoring
 # mem_monitor()

	# end loop over events
	Sndlhc.finish()
	
timer.Stop()
rtime = timer.RealTime()
ctime = timer.CpuTime()
print(' ') 
print("Real time ",rtime, " s, CPU time ",ctime,"s")
print(f'Digitised file written to {output_file}')
