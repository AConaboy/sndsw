#!/usr/bin/env python
import ROOT
import os,sys,subprocess
import time
import SndlhcGeo
from rootpyPickler import Unpickler

A,B=ROOT.TVector3(),ROOT.TVector3()

# for fixing a root bug,  will be solved in the forthcoming 6.26 release.
ROOT.gInterpreter.Declare("""
#include "MuFilterHit.h"
#include "AbsMeasurement.h"
#include "TrackPoint.h"

void fixRoot(MuFilterHit& aHit,std::vector<int>& key,std::vector<float>& value, bool mask) {
   std::map<int,float> m = aHit.GetAllSignals(false);
   std::map<int, float>::iterator it = m.begin();
   while (it != m.end())
    {
        key.push_back(it->first);
        value.push_back(it->second);
        it++;
    }
}
void fixRootT(MuFilterHit& aHit,std::vector<int>& key,std::vector<float>& value, bool mask) {
   std::map<int,float> m = aHit.GetAllTimes(false);
   std::map<int, float>::iterator it = m.begin();
   while (it != m.end())
    {
        key.push_back(it->first);
        value.push_back(it->second);
        it++;
    }
}
void fixRoot(MuFilterHit& aHit, std::vector<TString>& key,std::vector<float>& value, bool mask) {
   std::map<TString, float> m = aHit.SumOfSignals();
   std::map<TString, float>::iterator it = m.begin();
   while (it != m.end())
    {
        key.push_back(it->first);
        value.push_back(it->second);
        it++;
    }
}

void fixRoot(std::vector<genfit::TrackPoint*>& points, std::vector<int>& d,std::vector<int>& k, bool mask) {
      for(std::size_t i = 0; i < points.size(); ++i) {
        genfit::AbsMeasurement*  m = points[i]->getRawMeasurement();
        d.push_back( m->getDetId() );
        k.push_back( int(m->getHitId()/1000) );
    }
}
""")

Tkey  = ROOT.std.vector('TString')()
Ikey   = ROOT.std.vector('int')()
Value = ROOT.std.vector('float')()

class Monitoring():
   " set of monitor histograms "
   def __init__(self,options,FairTasks):

      self.simulation = options.simulation
      self.options = options
      self.EventNumber = -1
      self.TStart = -1
      self.TEnd   = -1
      self.Weight = 1

      path = options.path
      self.myclient = None

      if options.online:
         pass
      # setup geometry
      if self.simulation: self.snd_geo = SndlhcGeo.GeoInterface(options.path + options.geoFile)
      else: 
         if path.find('eos')>0: path  = options.server+options.path         
         self.snd_geo = SndlhcGeo.GeoInterface(path+options.geoFile)

      self.MuFilter = self.snd_geo.modules['MuFilter']
      self.Scifi    = self.snd_geo.modules['Scifi']
      self.systemAndPlanes = {1:2,2:5,3:7}
      self.zPos = self.getAverageZpositions()

      self.h = {}   # histogram storage

<<<<<<< HEAD
      self.runNr   = str(options.runNumber).zfill(6)
      self.FairTasks = {}
      for x in FairTasks:   #  keeps extended methods if from python class
         self.FairTasks[x.GetName()] = x
=======
        self.runNr   = str(options.runNumber).zfill(6)
<<<<<<< HEAD
# presenter file
<<<<<<< HEAD
        if hasattr(self, "saveTo") and options.saveTo!="":
          name = options.saveTo+'run'+self.runNr+'_'+str(options.nStart//1000000)+'.root'
        else:
           name = 'run'+self.runNr+'.root'
        if options.interactive: name = 'I-'+name
        self.presenterFile = ROOT.TFile(name,'recreate')
        for role in ['','shifter', 'expert']:
          self.presenterFile.mkdir('scifi/'+role)
          self.presenterFile.mkdir('mufilter/'+role)
          self.presenterFile.mkdir('daq/'+role)
          if options.interactive: self.presenterFile.mkdir('eventdisplay/'+role)
=======
      #   name = 'run'+self.runNr+'.root'
      #   if options.interactive: name = 'I-'+name
      #   self.presenterFile = ROOT.TFile(name,'recreate')
      #   self.presenterFile.mkdir('scifi')
      #   self.presenterFile.mkdir('mufilter')
      #   self.presenterFile.mkdir('daq')
      #   self.presenterFile.mkdir('eventdisplay')
>>>>>>> 3a2ccf43 (Updating files for use with simulation)
=======
>>>>>>> d8db0120 (at HEPHY)
        self.FairTasks = {}
        for x in FairTasks:   #  keeps extended methods if from python class
                 self.FairTasks[x.GetName()] = x
>>>>>>> cf0e3201a (Updating files for use with simulation)

<<<<<<< HEAD
<<<<<<< HEAD
      # setup input
      if self.simulation:
         partitions=[]

         eventChain = ROOT.TChain('cbmsim')

         if options.allFiles: 
            if options.simMode=='neutrino' and options.mode=='nue-extendedreconstruction':
               path=options.path+'sndlhc_13TeV_down_volTarget_100fb-1_SNDG18_02a_01_000/nueFilter/'
            else: 
               print(f'Options not set for all files and other modes')
               os.exit(1)
            allfiles=os.listdir(path)
            files= [ path+i for i in allfiles 
                     if all( [ i.split('.')[-1]=='root', i.find('geofile')==-1 ] ) ]
            [eventChain.Add(i) for i in files]

         elif not options.allFiles:
         # elif options.simMode in ('muonDIS', 'neutralhadron', 'neutrino', 'passingmuon'):
            files = [ options.path+options.fname ]
            eventChain.Add( files[0] ) # Only run when 1 file is used per job

         # else:
         #    allfiles = os.listdir(options.path)
         #    files = [ options.path+i for i in allfiles if all( [i.find(f"sndLHC.{options.simMode}")==0, i.find('digiCPP.root')>0] ) ]
         #    [eventChain.Add(f) for f in files]

      else:

         if options.customEventChain:
            eventChain=options.customEventChain
                     
            # Code added to analyse a collection of partitions from different runs e.g. for looking at mu_nu candidates 
            # Passing a dictionary of {runNr : partition}
            if options.signalpartitions:
               partitions=[f'sndsw_raw-{p}.root' for p in list(options.signalpartitions.values())]     

         else:
            partitions = []
            if path.find('eos')>0:
               # check for partitions
               dirlist  = str( subprocess.check_output("xrdfs "+options.server+" ls "+options.path+"run_"+self.runNr,shell=True) )
               for x in dirlist.split('\\n'):
                  ix = x.find('sndsw_raw-')
                  if ix<0: continue
                  partitions.append(x[ix:])
            else:
               # check for partitions
               dirlist  = os.listdir(options.path+"run_"+self.runNr)
               for x in dirlist:
                  if not x.find('sndsw_raw-')<0: partitions.append(x)
                  else: partitions = ["sndsw_raw-"+ str(options.partition).zfill(4)+".root"]

            if options.runNumber>0:
                  eventChain = ROOT.TChain('rawConv')
                  for p in partitions:
                     eventChain.Add(path+'run_'+self.runNr+'/'+p)

      rc = eventChain.GetEvent(0)
      if hasattr(eventChain, "EventHeader"):
         self.TStart = eventChain.EventHeader.GetEventTime()
         if options.nEvents <0:
            rc = eventChain.GetEvent(eventChain.GetEntries()-1)
         else:
            rc = eventChain.GetEvent(options.nEvents-1)
         self.TEnd = eventChain.EventHeader.GetEventTime()
=======
# setup input
        if options.online:
         pass 
        else:
            if self.simulation:
                eventChain = ROOT.TChain('cbmsim')
                allfiles = os.listdir(options.path)
                files = [ options.path+i for i in allfiles if all( [i.find(options.simEngine), i.find('dig.root')>0] ) ]
                [eventChain.Add(f) for f in files]
>>>>>>> d2ea8d70b (at HEPHY)
            
      # start FairRunAna
      self.run  = ROOT.FairRunAna()
      ioman = ROOT.FairRootManager.Instance()
      ioman.SetTreeName(eventChain.GetName())
      outFile = ROOT.TMemFile('dummy','CREATE')
      source = ROOT.FairFileSource(eventChain.GetCurrentFile()) # first file in chain is added here
      
      if self.simulation: 
         # Only need to add more files to source if multiple files are being used
         [source.AddFile(f) for f in files[1:]] # Skip first file

      elif options.customEventChain: # Code run when investigating numu candidates
         for idx, runNr in enumerate(options.signalpartitions):
            if idx!=0: source.AddFile(path+'run_'+runNr+'/'+partitions[idx]) # skip first partition which is added to the FairFileSource when it is instanced.

<<<<<<< HEAD
      else:
         for i in range(1,len(partitions)):
               p = partitions[i]
               source.AddFile(path+'run_'+self.runNr+'/'+p)
=======
            rc = eventChain.GetEvent(0)
            if hasattr(eventChain, "EventHeader"):
               self.TStart = eventChain.EventHeader.GetEventTime()
               if options.nEvents <0:
                  rc = eventChain.GetEvent(eventChain.GetEntries()-1)
               else:
                  rc = eventChain.GetEvent(options.nEvents-1)
               self.TEnd = eventChain.EventHeader.GetEventTime()
            
            # start FairRunAna
            self.run  = ROOT.FairRunAna()
            ioman = ROOT.FairRootManager.Instance()
            ioman.SetTreeName(eventChain.GetName())
            outFile = ROOT.TMemFile('dummy','CREATE')
            source = ROOT.FairFileSource(eventChain.GetCurrentFile())
            if self.simulation: 
               [source.AddFile(f) for f in files] 
            elif not options.customEventChain:
               for i in range(1,len(partitions)):
                     p = partitions[i]
                     source.AddFile(path+'run_'+self.runNr+'/'+p)
            else:
               for idx, runNr in enumerate(options.signalpartitions):
                  if idx!=0: source.AddFile(path+'run_'+runNr+'/'+partitions[idx]) # skip first partition which is added to the FairFileSource when it is instanced.
            self.run.SetSource(source)
            self.sink = ROOT.FairRootFileSink(outFile)
            self.run.SetSink(self.sink)
>>>>>>> d2ea8d70b (at HEPHY)

      self.run.SetSource(source)
      self.sink = ROOT.FairRootFileSink(outFile)
      self.run.SetSink(self.sink)

      for t in FairTasks: 
         self.run.AddTask(t)

      # avoiding some error messages
      xrdb = ROOT.FairRuntimeDb.instance()
      xrdb.getContainer("FairBaseParSet").setStatic()
      xrdb.getContainer("FairGeoParSet").setStatic()

<<<<<<< HEAD
      self.run.Init()
      if len(partitions)>0 or self.simulation:  self.eventTree = ioman.GetInChain()
      else:  self.eventTree = ioman.GetInTree()
=======
# initialize detector class for access to eventheader
        if not self.simulation:
         rc = eventChain.GetEvent(0)
         if not self.simulation:
            self.snd_geo.modules['Scifi'].InitEvent(eventChain.EventHeader)
            self.snd_geo.modules['MuFilter'].InitEvent(eventChain.EventHeader)
>>>>>>> d2ea8d70b (at HEPHY)
=======
      # setup input
      if self.simulation:
         partitions=[]
         
         #### For studying events that pass nue selection
         if options.fname == 'nue': options.fname='/eos/user/c/cvilela/SND_nue_analysis_May24/nuMC/filtered_stage1.root'

         eventChain = ROOT.TChain('cbmsim')
         if options.simTest or options.simMode in ('muonDIS', 'neutralhadron', 'neutrino', 'passingmuon', 'nue'):
            if options.simMode != 'nue': files = [ options.path+options.fname ]
            else: files = [ options.fname ]
            eventChain.Add( files[0] ) # Only run when 1 file is used per job
         else:
            allfiles = os.listdir(options.path)
            files = [ options.path+i for i in allfiles if all( [i.find(f"sndLHC.{options.simMode}")==0, i.find('digiCPP.root')>0] ) ]
            [eventChain.Add(f) for f in files]
            
      elif options.customEventChain:
         eventChain=options.customEventChain
                  
         # Code added to analyse a collection of partitions from different runs e.g. for looking at mu_nu candidates 
         # Passing a dictionary of {runNr : partition}
         if options.signalpartitions:
            partitions=[f'sndsw_raw-{p}.root' for p in list(options.signalpartitions.values())]

      else:
         partitions = []
         if path.find('eos')>0:
            # check for partitions
            dirlist  = str( subprocess.check_output("xrdfs "+options.server+" ls "+options.path+"run_"+self.runNr,shell=True) )
            for x in dirlist.split('\\n'):
               ix = x.find('sndsw_raw-')
               if ix<0: continue
               partitions.append(x[ix:])
         else:
            # check for partitions
            dirlist  = os.listdir(options.path+"run_"+self.runNr)
            for x in dirlist:
               if not x.find('sndsw_raw-')<0: partitions.append(x)
               else: partitions = ["sndsw_raw-"+ str(options.partition).zfill(4)+".root"]

         if options.runNumber>0:
               eventChain = ROOT.TChain('rawConv')
               for p in partitions:
                  eventChain.Add(path+'run_'+self.runNr+'/'+p)

      rc = eventChain.GetEvent(0)
      if hasattr(eventChain, "EventHeader"):
         self.TStart = eventChain.EventHeader.GetEventTime()
         if options.nEvents <0:
            rc = eventChain.GetEvent(eventChain.GetEntries()-1)
         else:
            rc = eventChain.GetEvent(options.nEvents-1)
         self.TEnd = eventChain.EventHeader.GetEventTime()
            
      # start FairRunAna
      self.run  = ROOT.FairRunAna()
      ioman = ROOT.FairRootManager.Instance()
      ioman.SetTreeName(eventChain.GetName())
      outFile = ROOT.TMemFile('dummy','CREATE')
      source = ROOT.FairFileSource(eventChain.GetCurrentFile()) # first file in chain is added here
      
      if self.simulation: 
         # Only need to add more files to source if multiple files are being used
         if not ( options.simTest or options.simMode in ('muonDIS', 'neutralhadron', 'neutrino', 'passingmuon', 'nue') ):
            [source.AddFile(f) for f in files[1:]] # Skip first file

      elif options.customEventChain: # Code run when investigating numu candidates
         for idx, runNr in enumerate(options.signalpartitions):
            if idx!=0: source.AddFile(path+'run_'+runNr+'/'+partitions[idx]) # skip first partition which is added to the FairFileSource when it is instanced.

      else:
         for i in range(1,len(partitions)):
               p = partitions[i]
               source.AddFile(path+'run_'+self.runNr+'/'+p)

      self.run.SetSource(source)
      self.sink = ROOT.FairRootFileSink(outFile)
      self.run.SetSink(self.sink)

      for t in FairTasks: 
         self.run.AddTask(t)

      # avoiding some error messages
      xrdb = ROOT.FairRuntimeDb.instance()
      xrdb.getContainer("FairBaseParSet").setStatic()
      xrdb.getContainer("FairGeoParSet").setStatic()

      self.run.Init()
      if len(partitions)>0 or self.simulation:  self.eventTree = ioman.GetInChain()
      else:  self.eventTree = ioman.GetInTree()

      # fitted tracks
      if "simpleTracking" in self.FairTasks:
         self.trackTask = self.FairTasks["simpleTracking"]
         self.Reco_MuonTracks = self.trackTask.fittedTracks
         self.clusMufi = self.trackTask.clusMufi
         self.clusScifi = self.trackTask.clusScifi
         self.trackTask.DSnPlanes = 3

      # initialize detector class for access to eventheader
      if not self.simulation:
         rc = eventChain.GetEvent(0)
         self.snd_geo.modules['Scifi'].InitEvent(eventChain.EventHeader)
         self.snd_geo.modules['MuFilter'].InitEvent(eventChain.EventHeader)
>>>>>>> dfd288d1f (loads of stuff on the simulation side and adding more details to scifi event t0)

<<<<<<< HEAD
      # fitted tracks
      if "simpleTracking" in self.FairTasks:
         self.trackTask = self.FairTasks["simpleTracking"]
         self.Reco_MuonTracks = self.trackTask.fittedTracks
         self.clusMufi = self.trackTask.clusMufi
         self.clusScifi = self.trackTask.clusScifi
         self.trackTask.DSnPlanes = 3

      # initialize detector class for access to eventheader
      if not self.simulation:
         rc = eventChain.GetEvent(0)
         self.snd_geo.modules['Scifi'].InitEvent(eventChain.EventHeader)
         self.snd_geo.modules['MuFilter'].InitEvent(eventChain.EventHeader)

         # get filling scheme, only necessary if not encoded in EventHeader, before 2022 reprocessing
         self.hasBunchInfo = False
         self.fsdict = False
         if hasattr(eventChain.EventHeader,"GetBunchType"):
            if not eventChain.EventHeader.GetBunchType()<0:
               self.hasBunchInfo = True
               print('take bunch info from event header')
         if not self.hasBunchInfo:
            try:
               fg  = ROOT.TFile.Open(options.server+options.path+'FSdict.root')
               pkl = Unpickler(fg)
               FSdict = pkl.load('FSdict')
               fg.Close()
               if options.runNumber in FSdict: self.fsdict = FSdict[options.runNumber]
            except:
               print('continue without knowing filling scheme',options.server+options.path)
=======
        # define the number of bunches in the LHC
        if eventChain.EventHeader.GetAccMode()==12: # ion runs
          self.Nbunches = 1782
          self.div = 8
        else: # proton runs
          self.Nbunches = 3564
          self.div = 4
        
        # get filling scheme, only necessary if not encoded in EventHeader, before 2022 reprocessing
        self.hasBunchInfo = False
        self.fsdict = False
        if hasattr(eventChain.EventHeader,"GetBunchType"):
           if not eventChain.EventHeader.GetBunchType()<0:
                self.hasBunchInfo = True
                print('take bunch info from event header')
        if not self.hasBunchInfo:
         try:
           fg  = ROOT.TFile.Open(options.server+options.path+'FSdict.root')
           pkl = Unpickler(fg)
           FSdict = pkl.load('FSdict')
           fg.Close()
           if options.runNumber in FSdict: self.fsdict = FSdict[options.runNumber]
         except:
           print('continue without knowing filling scheme',options.server+options.path)

         # get filling scheme, only necessary if not encoded in EventHeader, before 2022 reprocessing
         self.hasBunchInfo = False
         self.fsdict = False
         if hasattr(eventChain.EventHeader,"GetBunchType"):
            if not eventChain.EventHeader.GetBunchType()<0:
               self.hasBunchInfo = True
               print('take bunch info from event header')
         if not self.hasBunchInfo:
            try:
               fg  = ROOT.TFile.Open(options.server+options.path+'FSdict.root')
               pkl = Unpickler(fg)
               FSdict = pkl.load('FSdict')
               fg.Close()
               if options.runNumber in FSdict: self.fsdict = FSdict[options.runNumber]
            except:
               print('continue without knowing filling scheme',options.server+options.path)
<<<<<<< HEAD
>>>>>>> d8db0120 (at HEPHY)
         # if self.fsdict: 
         #   print('extract bunch info from filling scheme')
      #   if self.fsdict or self.hasBunchInfo: 
      #     for x in ['B1only','B2noB1','noBeam']:
      #        self.presenterFile.mkdir('mufilter/'+x)
      #        self.presenterFile.mkdir('scifi/'+x)
>>>>>>> 3a2ccf43 (Updating files for use with simulation)
>>>>>>> cf0e3201a (Updating files for use with simulation)
=======
>>>>>>> dfd288d1f (loads of stuff on the simulation side and adding more details to scifi event t0)

   def GetEntries(self):
       if  self.options.online:
         if  self.converter.newFormat:  return self.converter.fiN.Get('data').GetEntries()
         else:                                   return self.converter.fiN.event.GetEntries()
       else:
           return self.eventTree.GetEntries()

   def GetEvent(self,n):

      if self.eventTree.GetBranchStatus('Reco_MuonTracks'):
         for aTrack in self.eventTree.Reco_MuonTracks:
               if aTrack: aTrack.Delete()
         self.eventTree.Reco_MuonTracks.Delete()
      if "simpleTracking" in self.FairTasks:
         self.Reco_MuonTracks.Delete()
<<<<<<< HEAD

      self.eventTree.GetEvent(n)

      if not self.simulation:
         # initialize detector class for access to eventheader
         self.snd_geo.modules['Scifi'].InitEvent(self.eventTree.EventHeader)
         self.snd_geo.modules['MuFilter'].InitEvent(self.eventTree.EventHeader)

      if self.simulation: 
         self.Weight = self.eventTree.MCTrack[0].GetWeight()
      
      for t in self.FairTasks: 
            if t=='simpleTracking': self.FairTasks[t].ExecuteTask(nPlanes=3)
            else: self.FairTasks[t].ExecuteTask()

<<<<<<< HEAD
=======
      if self.options.online:
            online = self.options.online
            online.executeEvent(n)
            self.Reco_MuonTracks.Delete()
            if self.options.FairTask_convRaw:
                self.options.online.sTree.GetEvent(self.options.online.sTree.GetEntries()-1)
            for t in self.FairTasks: 
                if not t=='simpleTracking': self.FairTasks[t].ExecuteTask()
                else: self.FairTasks[t].ExecuteTask(options.trackType)
            self.eventTree = self.options.online.sTree
      else: 
            self.eventTree.GetEvent(n)
            if not self.simulation:
              # initialize detector class for access to eventheader
              self.snd_geo.modules['Scifi'].InitEvent(self.eventTree.EventHeader)
              self.snd_geo.modules['MuFilter'].InitEvent(self.eventTree.EventHeader)
            if self.simulation: self.Weight = self.eventTree.MCTrack[0].GetWeight()
            for t in self.FairTasks: 
                if t=='simpleTracking': self.FairTasks[t].ExecuteTask(nPlanes=3)
                else: self.FairTasks[t].ExecuteTask()
>>>>>>> d2ea8d70b (at HEPHY)
=======

      self.eventTree.GetEvent(n)

      if not self.simulation:
         # initialize detector class for access to eventheader
         self.snd_geo.modules['Scifi'].InitEvent(self.eventTree.EventHeader)
         self.snd_geo.modules['MuFilter'].InitEvent(self.eventTree.EventHeader)

      if self.simulation: 
         self.Weight = self.eventTree.MCTrack[0].GetWeight()
      
      for t in self.FairTasks: 
            if t=='simpleTracking': self.FairTasks[t].ExecuteTask(nPlanes=3)
            else: self.FairTasks[t].ExecuteTask()

>>>>>>> dfd288d1f (loads of stuff on the simulation side and adding more details to scifi event t0)
      self.EventNumber = n

# check for bunch xing type
      if not self.simulation:
         self.xing = {'all':True,'B1only':False,'B2noB1':False,'noBeam':False}
         if self.hasBunchInfo:
               binfo = self.eventTree.EventHeader
               self.xing['IP1']  = binfo.isIP1()
               self.xing['IP2']  = binfo.isIP2()
               self.xing['B1']   = binfo.isB1()
               self.xing['B2']   = binfo.isB2()
               self.xing['B1only']  = binfo.isB1Only()
               self.xing['B2noB1']  = binfo.isB2noB1()
               self.xing['noBeam']  = binfo.isNoBeam()
         elif self.fsdict:
               T   = self.eventTree.EventHeader.GetEventTime()
               bunchNumber = (T%(4*3564))//4
               nb1 = (3564 + bunchNumber - self.fsdict['phaseShift1'])%3564
               nb2 = (3564 + bunchNumber - self.fsdict['phaseShift1']- self.fsdict['phaseShift2'])%3564
               b1 = nb1 in self.fsdict['B1']
               b2 = nb2 in self.fsdict['B2']
               IP1 = False
               IP2 = False
               if b1:
                  IP1 =  self.fsdict['B1'][nb1]['IP1']
               if b2:
                  IP2 =  self.fsdict['B2'][nb2]['IP2']
               self.xing['IP1']  = IP1
               self.xing['IP2']  = IP2
               self.xing['B1']   = b1
               self.xing['B2']   = b2
               self.xing['B1only']   = b1 and not IP1 and not b2
               self.xing['B2noB1']  = b2 and not b1
               self.xing['noBeam'] = not b1 and not b2
               if self.xing['B1only']  and self.xing['B2noB1']  or self.xing['B1only'] and self.xing['noBeam'] : print('error with b1only assignment',self.xing)
               if self.xing['B2noB1']  and self.xing['noBeam'] : print('error with b2nob1 assignment',self.xing)

      return self.eventTree

<<<<<<< HEAD
<<<<<<< HEAD
=======
<<<<<<< HEAD
   def publishRootFile(self):
   # try to copy root file with TCanvas to EOS
       self.presenterFile.Close()
       if self.options.online:
           wwwPath = "/eos/experiment/sndlhc/www/online"
       elif self.options.path.find('2022')>0 :
           wwwPath = "/eos/experiment/sndlhc/www/reprocessing"
       else:    
           wwwPath = "/eos/experiment/sndlhc/www/offline"
       if self.options.sudo: 
         try:
            rc = os.system("xrdcp -f "+self.presenterFile.GetName()+"  "+os.environ['EOSSHIP']+wwwPath)
         except:
            print("copy of root file failed. Token expired?")
         self.presenterFile = ROOT.TFile('run'+self.runNr+'.root','update')

   def purgeMonitorHistos(self):
        wwwPath = "/eos/experiment/sndlhc/www/online/"
        for r in self.options.runNumbers.split(','):
            if r!= '': runList.append(int(r))
        for r in runList:
              runNr   = str(r).zfill(6)+'.root'
              f = ROOT.TFile.Open(options.server+wwwPath+ runNr)
              g = ROOT.TFile('tmp.root','recreate')
              for d in f.GetListOfKeys():
                  D = f.Get(d.GetName())
                  D.Purge()
              for d in f.GetListOfKeys():
                name = d.GetName()
                s = f.Get(name)
                g.mkdir(name)
                g.cd(name)
                for x in s.GetListOfKeys():
                    s.Get(x.GetName()).Write()
        g.Close()
        f.Close()
        os.system('xrdcp -f tmp.root  '+self.options.server+options.path+rName)

   def updateHtml(self):
      if self.options.online: destination="online"
      elif self.options.path.find('2022')>0: destination="reprocessing"
      else: destination="offline"
      rc = os.system("xrdcp -f "+os.environ['EOSSHIP']+"/eos/experiment/sndlhc/www/"+destination+".html  . ")
      old = open(destination+".html")
      oldL = old.readlines()
      old.close()
      tmp = open("tmp.html",'w')
      found = False
      for L in oldL:
           if not L.find(self.runNr)<0: return
           if L.find("https://snd-lhc-monitoring.web.cern.ch/"+destination+"/run.html?file=run")>0 and not found:
              r = str(self.options.runNumber)
              Lnew = '            <li> <a href="https://snd-lhc-monitoring.web.cern.ch/'+destination+'/run.html?file=run'
              Lnew+= self.runNr+'.root&lastcycle">run '+r+' </a>'+self.options.startTime
              if self.options.postScale>1: Lnew+="  post scaled:"+str(self.options.postScale)
              Lnew+='\n'
              tmp.write(Lnew)
              found = True
           tmp.write(L)
      tmp.close()
      os.system('cp '+destination+'.html '+destination+time.ctime().replace(' ','')+'.html ')  # make backup
      os.system('cp tmp.html '+destination+'.html')
      try:
            rc = os.system("xrdcp -f "+destination+".html  "+os.environ['EOSSHIP']+"/eos/experiment/sndlhc/www/")
      except:
            print("copy of html failed. Token expired?")

   def cleanUpHtml(self,destination="online"):
      rc = os.system("xrdcp -f "+os.environ['EOSSHIP']+"/eos/experiment/sndlhc/www/"+destination+".html  . ")
      old = open(destination+".html")
      oldL = old.readlines()
      old.close()
      tmp = open("tmp.html",'w')
      dirlist  = str( subprocess.check_output("xrdfs "+os.environ['EOSSHIP']+" ls /eos/experiment/sndlhc/www/"+destination+"/",shell=True) ) 
      for L in oldL:
           OK = True
           if L.find("https://snd-lhc-monitoring.web.cern.ch/"+destination)>0:
              k = L.find("file=")+5
              m =  L.find(".root")+5
              R = L[k:m]
              if not R in dirlist: OK = False  
           if OK: tmp.write(L)
      tmp.close()
      os.system('cp tmp.html '+destination+'.html')
      try:
            rc = os.system("xrdcp -f "+destination+".html  "+os.environ['EOSSHIP']+"/eos/experiment/sndlhc/www/")
      except:
            print("copy of html failed. Token expired?")

   def checkAlarm(self,minEntries=1000):
      self.alarms = {'scifi':[]}
      # check for empty histograms in scifi signal
      entries = {}
      for mat in range(30):
           entries[mat] = self.h['scifi-mat_'+str(mat)].GetEntries()
      res = sorted(entries.items(), key=operator.itemgetter(1),reverse=True)
      if res[0][1]>minEntries: # choose limit which makes it sensitive to check for empty mats
          for mat in range(30):
              if entries[mat] < 1:
                   s = mat//6 + 1
                   if mat%6 < 3 : p='H'
                   else : p='V'
                   e = str(s)+p+str(mat%3)
                   self.alarms['scifi'].append(e)
      if len(self.alarms['scifi']) > 0:  self.sendAlarm()

      self.alarms['Veto']=[]
      # check for empty histograms in veto signal
      entries = {}
      s = 1
      for p in range(2):
         for side in ['L','R']:
             entries[str(10*s+p)+side] = self.h['mufi-sig'+p+str(10*s+p)].GetEntries()
      res = sorted(entries.items(), key=operator.itemgetter(1),reverse=True)
      if res[0][1]>minEntries: # choose limit which makes it sensitive to check for empty boards
       for p in range(2):
         for side in ['L','R']:
             if entries[str(10*s+p)+side] < 1:
                self.alarms['Veto'].append(str(10*s+p)+side())
      # check for empty histograms in US signal
      entries = {}
      s = 2
      for p in range(5):
         for side in ['L','R']:
             entries[str(s*10+p)+side] = self.h['mufi-sig2'+p+str(s*10+p)].GetEntries()
      res = sorted(entries.items(), key=operator.itemgetter(1),reverse=True)
      if res[0][1]>1000: # choose limit which makes it sensitive to check for empty boards
       for p in range(5):
         for side in ['L','R']:
             if entries[str(s*10+p)+side] < 1:
                self.alarms['Veto'].append(str(s*10+p)+side())


   def sendAlarm(self):
       print('ALARM',self.alarms)

=======
>>>>>>> d8db0120 (at HEPHY)
>>>>>>> d2ea8d70b (at HEPHY)
=======
>>>>>>> dfd288d1f (loads of stuff on the simulation side and adding more details to scifi event t0)
   def systemAndOrientation(self,s,plane):
      if s==1 or s==2: return "horizontal"
      if plane%2==1 or plane == 6: return "vertical"
      return "horizontal"

   def getAverageZpositions(self):
      zPos={'MuFilter':{},'Scifi':{}}
      for s in self.systemAndPlanes:
          for plane in range(self.systemAndPlanes[s]):
             bar = 4
             p = plane
             if s==3 and (plane%2==0 or plane==7): 
                bar = 90
                p = plane//2
             elif s==3 and plane%2==1:
                bar = 30
                p = plane//2
             self.MuFilter.GetPosition(s*10000+p*1000+bar,A,B)
             zPos['MuFilter'][s*10+plane] = (A.Z()+B.Z())/2.
      for s in range(1,6):
         mat   = 1
         sipm = 1
         channel = 64
         for o in range(2):
             self.Scifi.GetPosition(channel+1000*sipm+10000*mat+100000*o+1000000*s,A,B)
             zPos['Scifi'][s*10+o] = (A.Z()+B.Z())/2.
      return zPos

#  Scifi specific code
   def Scifi_xPos(self,detID):
        orientation = (detID//100000)%10
        nStation = 2*(detID//1000000-1)+orientation
        mat   = (detID%100000)//10000
        X = detID%1000+(detID%10000)//1000*128
        return [nStation,mat,X]   # even numbers are Y (horizontal plane), odd numbers X (vertical plane)
