#!/usr/bin/env python
import ROOT, os, csv, SndlhcGeo

class Numusignaleventtiming(object):

    # def __init__(self, options, tw):
    def Init(self, options, monitor):
               
        self.options=options
        self.tw=tw
        self.M=tw.M
        self.runNr = tw.runNr 
        self.muAna = tw.muAna
        self.timealignment=tw.timealignment

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

        self.freq=160.316E6
        self.TDC2ns=1E9/self.freq

        self.xrefs=tw.xrefs
        self.sides=('left', 'right')

        self.hists=tw.hists
        
        self.sigmatds0=0.263 # ns 

        if not hasattr(self.muAna, 'task'): self.muAna.SetTask(tw)

    def InvestigateSignalEvents(self):

        numusignalevent_filepath = '/afs/cern.ch/work/a/aconsnd/numusignalevents.csv'
        with open(numusignalevent_filepath, 'r') as f:
            reader=csv.reader(f)
            nu_mu_data=[r for r in reader]
        self.nu_mu_events={int(x[0]):int(x[1]) for x in nu_mu_data}

        self.eventChains={runNr:ROOT.TChain('rawConv') for runNr in self.nu_mu_events} 
        
        for runNr in self.nu_mu_events:
            eventNumber=self.nu_mu_events[runNr]

            runNumber = str(runNr).zfill(6)
            partition = int(eventNumber) // int(1E6)
            dirlist  = os.listdir(self.options.path+"run_"+runNumber)
            partitionfile=f'{self.options.path}run_{runNumber}/sndsw_raw-{str(partition).zfill(4)}.root'
            self.eventChains[runNr].Add(partitionfile)

        self.SignalEventTimingHists()  

    def SignalEventTimingHists(self):

        for runNr in self.eventChains:
            
            chain = self.eventChains[runNr]
            chain.GetEvent(0)

            self.muAna.GetRunYear(runNr)
            geofile=self.muAna.GetGeoFile(runNr)
            self.M.eventTree = chain
            self.M.snd_geo = SndlhcGeo.GeoInterface(self.options.server + self.options.path + geofile)
            self.M.MuFilter = self.M.snd_geo.modules['MuFilter']
            self.M.Scifi = self.M.snd_geo.modules['Scifi']
            self.M.zPos = self.M.getAverageZpositions()

            # Load signal event
            signalevent=self.nu_mu_events[runNr]
            chain.GetEvent(signalevent)

            ### Evaluate tracking
            self.M.FairTasks["simpleTracking"].ExecuteTask(nPlanes=3)

            tracks={1:[], 3:[]}
            Reco_MuonTracks=self.M.Reco_MuonTracks
            inVeto, inDS=False, False
            for i,track in enumerate(Reco_MuonTracks):
                if any([not track.getFitStatus().isFitConverged(), track.getFitStatus().getNdf()==0]): continue
                if track.GetUniqueID()==1: 
                    inVeto=True
                    tracks[1].append(Reco_MuonTracks[i])
                if track.GetUniqueID()==3: 
                    inDS=True
                    tracks[3].append(Reco_MuonTracks[i])
            if not inDS: return

            ### If there are more than 1 DS track, take the track with the lowest chi2/Ndf
            if len(tracks[3])==1: dstrack= tracks[3][0]
            else: 
                tmp={i:i.getFitStatus().getChi2()/i.getFitStatus().getNdf() for i in tracks[3]}
                dstrack=tmp[ min(tmp) ]
            
            fitStatus=dstrack.getFitStatus()
            fstate=dstrack.getFittedState()

            self.tw.pos=fstate.getPos()
            self.tw.mom=fstate.getMom()
            self.trackchi2NDF=fitStatus.getChi2()/fitStatus.getNdf()+1E-10            

            # Get DS event t0
            x=self.muAna.GetDSHaverage(event.Digi_MuFilterHits) # Now returns in ns
            if x==-999.: 
                print(f'Event {self.M.EventNumber} has no DS horizontal hits with 2 fired SiPMs')
                return
            self.TDS0, firedDSHbars=x

            hits=chain.Digi_MuFilterHits
            self.tw.sp.FillHists(hits)