#!/usr/bin/env python
import ROOT, csv, os, pickle
import numpy as np
from pathlib import Path
from HCALTools import HCALTools 
ROOT.gInterpreter.ProcessLine('#include "/afs/cern.ch/user/a/aconsnd/sndsw/analysis/tools/sndSciFiTools.h"')

class ExtendedMuonReconstruction(object):

    def __init__(self, options, tw):
        self.options=options
        self.tw=tw
        self.simulation=tw.simulation
        self.runNr = tw.runNr 
        self.muAna = tw.muAna
        self.timealignment=tw.timealignment
        self.TWCorrectionRun=tw.TWCorrectionRun

        ### If no time-walk correction run is provided. Set the default correction run depending on time alignment of the data set
        self.filekey = options.fname.split('/')[-1].split('_')[-1].replace('.root', '')

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
        self.Scifi = tw.Scifi
        self.barlengths = self.muAna.BuildBarLengths(self.MuFilter)

        self.freq=160.316E6
        self.TDC2ns=1E9/self.freq

        self.sides=('left', 'right')

        self.hists=tw.hists

        self.sigmatds0=0.263 # ns 

        self.A, self.B = ROOT.TVector3(), ROOT.TVector3()
        self.MuFilter.GetPosition(24004, self.B, self.B)
        self.HCAL5z = 0.5*(self.A.z() + self.B.z())

        self.hcalTools = HCALTools(self.muAna, self.MuFilter)
        self.hcalTools.filekey=self.filekey

        """
        Set xy acceptance limits for last target wall
        this will help balance the dataset passed to the BDT.
        Not hugely urgent but I should do it.
        """
        # dummy_detID=int(1e6*5 + 1e5*1 + 1e4 + 1e3 + 1)
        # self.Scifi.GetSiPMPosition(, self.B, self.B)
        # self.HCAL5z = 0.5*(self.A.z() + self.B.z())        
        
        if options.signalpartitions: self.Loadnumuevents()

        self.numuStudy=True if options.numuStudy else False 

        # Not exact, just rough for rejecting rubbish combinations of DS clusters
        self.acceptancelimits={'x':[-100, 20], 'y':[-10, 100]}

        if self.options.mode=='nue-extendedreconstruction' and self.simulation: 

            d = f'{self.outpath}{self.tw.mode}/'
            os.makedirs(d, exist_ok=True)

            dirkey1, dirkey2, filename = self.options.fname.split('/')
            key=filename.replace('.root', '').split('_')[1]

            self.keynamedict = {'wMuon':'with muon', 'woMuon': 'w/o muon', 'allEvents':'all events'}

            self.hcalTools.datafilename=d+f'extendedreconstruction_{key}.csv'
            column_names=['filekey', 'EventNumber','hasMuon','interactionWall',
            'scifi_median_x','scifi_median_y',
            'scifi_residual_x','scifi_residual_y',
            'dx0','dx1','dx2','dx3','dx4',
            'dy0','dy1','dy2','dy3','dy4',
            'ds0','ds1','ds2','ds3','ds4', 'ds_scifi',
            'x0','x1','x2','x3','x4',
            'y0','y1','y2','y3','y4',
            'lambdax0','lambdax1', 'lambdax2','lambdax3', 'lambdax4', 
            'lambday0','lambday1', 'lambday2','lambday3', 'lambday4',
            'HCAL5barcode'
            ]

            with open(self.hcalTools.datafilename, 'w') as f:
                writer=csv.writer(f)
                writer.writerow(column_names)

        elif self.options.mode=='nue-extendedreconstruction' and not self.simulation: 

            # d = f'{self.outpath}{self.tw.mode}/'
            # os.makedirs(d, exist_ok=True)

            # dirkey1, dirkey2, filename = self.options.fname.split('/')
            # key=filename.replace('.root', '').split('_')[1]

            self.keynamedict = {'wMuon':'with muon', 'woMuon': 'w/o muon', 'allEvents':'all events'}

            # self.datafilename=d+f'extendedreconstruction_{key}.csv'
            self.column_names=['filekey', 'EventNumber', 'hasMuon', 'interactionWall',
            'scifi_median_x','scifi_median_y',
            'scifi_residual_x','scifi_residual_y',
            'dx0','dx1','dx2','dx3','dx4',
            'dy0','dy1','dy2','dy3','dy4',
            'x0','x1','x2','x3','x4',
            'y0','y1','y2','y3','y4',
            'lambdax0','lambdax1', 'lambdax2','lambdax3', 'lambdax4', 
            'lambday0','lambday1', 'lambday2','lambday3', 'lambday4',
            'HCAL5bars'
            ]

    def ExtendReconstruction(self, hits, scifi_hits, mode='write'):
        # Here I want to get the points in space from the DS hits, and see if a US hit aligns with these
        # If they do then I can plot the doca between the line formed between these DS hits and the US hit
        # if self.options.OutgoingMuon=='yes' and not eventHasMuon: return
        # elif self.options.OutgoingMuon=='no' and eventHasMuon: return
        if self.simulation:
            self.hcalTools.eventHasMuon=self.hcalTools.OutgoingMuon(self.tw.M.eventTree)
            if self.hcalTools.eventHasMuon: self.muonhistkey = 'wMuon'
            elif not self.hcalTools.eventHasMuon: self.muonhistkey = 'woMuon'
        self.hcalTools.EventNumber = self.tw.M.EventNumber

        self.hcalTools.barycentres = self.muAna.GetBarycentres(hits, MuFilter=self.MuFilter)
        self.hcalTools.xbarycentres = self.muAna.GetOverallXBarycentre(self.hcalTools.barycentres, mode='maxQDC')

        self.hcalTools.dsCluster(hits, self.MuFilter)
        if len(self.hcalTools.clusMufi)==0: 
            return

        self.hcalTools.GetDSClusterCentroids()
        
        fired_planes=list(self.hcalTools.DS_centroids.keys())
        if len(fired_planes)==3: 
            print(f'3 fired DS planes in event {self.tw.M.EventNumber}')
            return     
        if len(fired_planes)==0:
            return 

        histname = 'n_DSclusters'
        if not histname in self.hists:
            title="Number of cluster formed in muon system;# DS clusters;Counts"
            self.hists[histname] = ROOT.TH1F(histname, title, 26, 0, 26)
        self.hists[histname].Fill(len(self.hcalTools.DS_centroids))
        
        """
        For each permutation of pairs of x,y fired bars, I can make a straight line
        and extrapolate to the US to check the doca with the x-barycentre, y-barycentre.

        In order to know which permutation is successful, I will need to keep track of the detector IDs somehow. 
        """
        
        self.hcalTools.GetCombinatorics()

        for p in ('x', 'y'):
            histname=f'n-{p}z-combinations'
            if not histname in self.hists:
                title=f'Number of {p} clusters permutations in the DS;# {p} cluster permutations;Counts'
                self.hists[histname] = ROOT.TH1F(histname, title, 6, 0, 6)
            self.hists[histname].Fill(len(self.hcalTools.combinations[p]))

        # Count HCAL hits in each plane
        # if mode=='investigate': self.GetMultiplicity(hits)

        self.hcalTools.interactionWall = self.hcalTools.GetInteractionWall(scifi_hits)
        self.hcalTools.Get_interaction_median_positions(scifi_hits, self.Scifi)

        # If 2 fired planes in the DS, I can connect the space points in each combination together and look for a hit in the US
        
        if len(fired_planes)==2:
            self.hcalTools.xy_residuals = {plane:{'x':np.nan, 'y':np.nan} for plane in range(5)}
            self.hcalTools.xy_residuals['intWall']={'x':np.nan, 'y':np.nan}

            for idx, proj in enumerate(['x', 'y']):
                
                residuals=[]
                for combination in self.hcalTools.combinations[proj]:
                    
                    # Make lines for the xz and yz projections that join the points of this pair
                    line = self.hcalTools.ConnectPoints(combination, proj)
                    # line == False if the points draw a line out of the acceptance of HCAL plane 5
                    if not line: continue

                    # returns a dictionary of the residual in that projection
                    res = self.hcalTools.USresidual(line, proj,self.MuFilter,self.Scifi) 
                    if not 4 in res and 5 in res: continue

                    residuals.append(res)

                if len(residuals)==0: continue

                # Define best combination in each projection as the one with the lowest residual. Projections are orthogonal so no issue.
                best_residual  = min(residuals, key=lambda d:d['intWall'])

                for plane in best_residual:
                    # if plane=='intWall':continue
                    # Update value for each projection if there are suitable combinations
                    self.hcalTools.xy_residuals[plane][proj] = best_residual[plane]
 
            # Require that the xy_residual is defined for the 4th and 5th plane
            if list(self.hcalTools.xy_residuals[4].values()) == [np.nan, np.nan]: 
                print(f'Event {self.tw.M.EventNumber}, no x and y residual in plane 5')
                return 
            if list(self.hcalTools.xy_residuals[3].values()) == [np.nan, np.nan]: 
                print(f'Event {self.tw.M.EventNumber}, no x and y residual in plane 4')
                return

            self.hcalTools.lambda_x_dict = {i:np.nan for i in range(5)}
            self.hcalTools.lambda_y_dict = {i:np.nan for i in range(5)}

            for plane in self.hcalTools.xy_residuals:
                if plane=='intWall':continue
                lambda_x = abs(self.hcalTools.xbarycentres[plane]['lambda_x'])
                self.hcalTools.lambda_x_dict[plane] = lambda_x
                lambda_y = self.hcalTools.barycentres[plane]['y-barycentre']['lambda_y']
                self.hcalTools.lambda_y_dict[plane] = lambda_y

                if len(self.hcalTools.xy_residuals[plane])==2:

                    self.xyresiduals_hists(plane)
                    # self.USmultvds_hists(plane)

                if mode=='write': 
                    self.lambda_hists(plane)
                    self.ds_hists(plane)

            self.hcalTools.GetHCAL5barscode(hits)
            self.hcalTools.Get_ds()

            self.hcalTools.getdata(mode=mode)

        elif len(fired_planes)==1:
            pass
        else: pass

    def GetHCAL5barscode(self, hits):
        c=[False]*10 
        US_detID_list = [i.GetDetectorID() for i in hits if all([self.muAna.parseDetID(i.GetDetectorID())[0]==2, self.muAna.parseDetID(i.GetDetectorID())[1]==4])]
        for detID in US_detID_list:
            bar=self.muAna.parseDetID(detID)[2]
            c[bar]=True 
        
        self.HCAL5barscode=''.join(['1' if state else '0' for state in c])

    def OutgoingMuon(self):
        ### Return true for events with an outgoing muon
        if not self.simulation:
            print(f'No MCTrack branch in real data! ')
            return False
        if abs(self.tw.M.eventTree.MCTrack[1].GetPdgCode()) == 13: return True 
        else: return False

    def xyresiduals_hists(self, plane):

        for key in ('wMuon', 'woMuon', 'allEvents'):
            histname=f'xyresidual_{key}_plane{plane}'
            if not histname in self.hists:
                title='#splitline{Correlation of residuals '+self.keynamedict[key]+'}{between DS hits and barycentre in plane '+str(plane+1)+'};#Delta(DS hits, x-barycentre) [cm];#Delta(DS hits, y-barycentre) [cm];Counts'
                self.hists[histname]=ROOT.TH2F(histname, title, 100, -50, 50, 100, -50, 50)
        histname=f'xyresidual_{self.muonhistkey}_plane{plane}'
        self.hists[histname].Fill(*self.hcalTools.xy_residuals[plane].values())                
        self.hists[f'xyresidual_allEvents_plane{plane}'].Fill(*self.hcalTools.xy_residuals[plane].values())

    def lambda_hists(self, plane):

        for key in ('wMuon', 'woMuon', 'allEvents'):

            histname=f'lambdax_{key}_plane{plane}'
            if not histname in self.hists:
                title='#splitline{Energy deposition width in x in plane '+str(plane+1)+'}{'+self.keynamedict[key]+'};#lambda_{x} = x_{L} - x_{R} [cm];Counts'
                self.hists[histname]=ROOT.TH1F(histname, title, 50, 0, 50)

            histname=f'lambday_{key}_plane{plane}'
            if not histname in self.hists:
                title='#splitline{Energy deposition width in y in plane '+str(plane+1)+'}{'+self.keynamedict[key]+'};#lambda_{y} = #sum_{bars} #frac{QDC_{bar}}{QDC_{plane}}y_{bar} [cm];Counts'
                self.hists[histname]=ROOT.TH1F(histname, title, 50, 0, 50)

        histname = f'lambdax_{self.muonhistkey}_plane{plane}'
        lambda_x = self.hcalTools.lambda_x_dict[plane]
        if lambda_x: 
            self.hists[histname].Fill(lambda_x)
            self.hists[f'lambdax_allEvents_plane{plane}'].Fill(lambda_x)
        
        histname = f'lambday_{self.muonhistkey}_plane{plane}'
        lambda_y = self.hcalTools.lambda_y_dict[plane]
        if lambda_y: 
            self.hists[histname].Fill(lambda_y)
            self.hists[f'lambday_allEvents_plane{plane}'].Fill(lambda_y)

    def ds_hists(self, plane):

        for key in ('wMuon', 'woMuon', 'allEvents'):

            histname=f'dx_{key}_plane{plane}'
            if not histname in self.hists:
                title='#splitline{x residual between DS hits and expected position in plane '+str(plane+1)+'}{'+self.keynamedict[key]+'};dx^{2} [cm];Counts'
                self.hists[histname]=ROOT.TH1F(histname, title, 50, 0, 50)

            histname=f'dy_{key}_plane{plane}'
            if not histname in self.hists:
                title='#splitline{y residual between DS hits and expected position in plane '+str(plane+1)+'}{'+self.keynamedict[key]+'};dy^{2} [cm];Counts'
                self.hists[histname]=ROOT.TH1F(histname, title, 50, 0, 50)

            histname=f'ds_{key}_plane{plane}'
            if not histname in self.hists:
                title='Resultant residual between DS hits and expected position in plane '+str(plane+1)+'}{'+self.keynamedict[key]+'};ds = #sqrt{dx^{2} + dy^{2}} [cm];Counts'
                self.hists[histname]=ROOT.TH1F(histname, title, 50, 0, 50)
        
        if 'x' in self.hcalTools.xy_residuals[plane]:
            histname = f'dx_{self.muonhistkey}_plane{plane}'
            self.hists[histname].Fill(self.hcalTools.xy_residuals[plane]['x'])
            self.hists[f'dx_allEvents_plane{plane}'].Fill(self.hcalTools.xy_residuals[plane]['x'])
        
        if 'y' in self.hcalTools.xy_residuals[plane]:
            histname = f'dy_{self.muonhistkey}_plane{plane}'
            self.hists[histname].Fill(self.hcalTools.xy_residuals[plane]['y'])
            self.hists[f'dy_allEvents_plane{plane}'].Fill(self.hcalTools.xy_residuals[plane]['y'])
        
        if 'x' in self.hcalTools.xy_residuals[plane] and 'y' in self.hcalTools.xy_residuals[plane]:
            dx, dy = self.hcalTools.xy_residuals[plane]['x'], self.hcalTools.xy_residuals[plane]['y']
            ds = np.sqrt(dx**2 + dy**2)
            histname = f'ds_{self.muonhistkey}_plane{plane}'
            self.hists[histname].Fill(ds)
            self.hists[f'ds_allEvents_plane{plane}'].Fill(ds)

    def GetLambda(self, plane, proj):
        if proj=='x': 
            if plane not in self.hcalTools.xbarycentres: return
            if not 'lambda_x' in self.hcalTools.xbarycentres[plane]:return
            b=self.hcalTools.xbarycentres[plane]['lambda_x']
        elif proj=='y': 
            if plane not in self.hcalTools.barycentres: return
            if not 'y-barycentre' in self.hcalTools.barycentres[plane]:return
            if not 'lambda_y' in self.hcalTools.barycentres[plane]['y-barycentre']:return
            b=self.hcalTools.barycentres[plane]['y-barycentre']['lambda_y']
        return b        

    def GetMultiplicity(self, hits):
        
        self.multiplicity_dict = {2:{i:0 for i in range(5)}, 3:{i:0 for i in range(7)}}
        
        for hit in hits:
            detID=hit.GetDetectorID()
            s,p,b=self.muAna.parseDetID(detID)

            if s==2: self.multiplicity_dict[s][p]+=1
            if s==3: 
                DSplanenumber=self.muAna.GetDSPlaneNumber(detID)
                self.multiplicity_dict[s][DSplanenumber]+=1

        for key in ('wMuon', 'woMuon', 'allEvents'):

            for plane in range(5):
                histname=f'USmultiplicity_{key}_plane{plane}'
                if not histname in self.hists:
                    title='#splitline{Number of fired bars in HCAL plane '+str(plane+1)+'}{'+self.keynamedict[key]+'};N fired bars;Counts'
                    self.hists[histname]=ROOT.TH1I(histname, title, 11, 0, 11)
        if self.simulation:                
            if self.eventHasMuon: self.hists[f'USmultiplicity_wMuon_plane{plane}'].Fill(self.multiplicity_dict[2][plane])
            elif not self.eventHasMuon: self.hists[f'USmultiplicity_woMuon_plane{plane}'].Fill(self.multiplicity_dict[2][plane])
        self.hists[f'USmultiplicity_allEvents_plane{plane}'].Fill(self.multiplicity_dict[2][plane])

    def WriteOutHistograms(self):

        if not self.simulation:
            d = f'{self.outpath}splitfiles/run{self.runNr}/{self.tw.mode}/'
            os.makedirs(d, exist_ok=True)
            
            outfilename=d+f'extendedreconstruction_{self.options.nStart}.root' 

        elif self.simulation and self.options.simMode=='neutrino':

            if self.options.mode=='nue-extendedreconstruction': 

                d = f'{self.outpath}{self.tw.mode}/'
                os.makedirs(d, exist_ok=True)

                dirkey1, dirkey2, filename = self.options.fname.split('/')
                key=filename.replace('.root', '').split('_')[1]

                # if self.options.OutgoingMuon=='yes':   muonkey='wMuon'
                # elif self.options.OutgoingMuon=='no':   muonkey='woMuon'
                # elif self.options.OutgoingMuon=='all':   muonkey='allEvents'

                outfilename = d+f'extendedreconstruction_{key}.root'

            # d = f'{self.outpath}{self.tw.mode}/'
            # os.makedirs(d, exist_ok=True)

            # dirkey1, dirkey2, filename = self.options.fname.split('/')
            # key=filename.split('_')[1]
            # outfilename=d+f'extendedreconstruction_{key}.root'

        elif self.simulation and self.options.simMode == 'neutralhadron':
            d = f'{self.outpath}{self.tw.mode}/'
            os.makedirs(d, exist_ok=True)

            particle_type, Emin, Emax, key = self.options.fname.split('_')[3:7]
            outfilename=d+f'extendedreconstruction_{particle_type}_{Emin}_{Emax}_{key}.root'

        elif self.simulation and self.options.simMode == 'passingmuon':
            d = f'{self.outpath}{self.tw.mode}/'
            os.makedirs(d, exist_ok=True)

            keys=self.options.fname.split('_')[1:3]
            outfilename=d+f'extendedreconstruction_{keys[0]}_{keys[1]}.root'

        elif self.simulation and self.options.simMode == 'nue':
            print(f'Not implemented nue saving protocol! ')     

        if os.path.exists(outfilename): outfile=ROOT.TFile.Open(outfilename, 'recreate')
        else: outfile=ROOT.TFile.Open(outfilename, 'create')            

        for hname in self.hists:
            
            hist=self.hists[hname]
            if hname in ('yEx', 'xEx', 'reft', 'n-xz-combinations', 'n-yz-combinations', 'n_DSclusters', 'USDSmultiplicityvresresidual'):
                outfile.WriteObject(hist, hname, 'kOverwrite')
                continue

            key, muonkey, plane = hname.split('_')

            if not hasattr(outfile, key): folder=outfile.mkdir(key)
            else: folder=outfile.Get(key)
            
            folder.cd()
            hist.Write(hname, 2) # The 2 means it will overwrite a hist of the same name            

        outfile.Close()
        print(f'{len(self.hists)} histograms saved to {outfilename}')   

        print(f'Data written to {self.hcalTools.datafilename}') 

class EMRresults(object):
    def __init__(self):
        self.outpath = '/eos/experiment/sndlhc/users/aconsnd/simulation/neutrino/data/nue-extendedreconstruction/'
        self.muAna = SetUpAnalysisClass()
        self.flags=['woMuon', 'wMuon', 'allEvents']

    def GetHists(self):

        files = {flag:f'{self.outpath}extendedreconstruction_{flag}.root' for flag in self.flags}

        self.hists={flag:{} for flag in self.flags}
        for flag in files: 
            f=ROOT.TFile.Open(files[flag])
            for plane in range(5):
                hist=f.Get(f'resultantresidual_plane{plane}')
                hist.SetDirectory(ROOT.gROOT)
                self.hists[flag][hist.GetName()]=hist
            f.Close()

    def OverlayHists(self, modes=['allEvents', 'wMuon', 'woMuon'], histname='resultantresidual_plane4'):
        c=ROOT.TCanvas()
        c.SetTitle(f'Comparison of histograms')

        l = ROOT.TLegend()
        colours=[ROOT.kRed, ROOT.kBlack, ROOT.kBlue]

        for i, flag in enumerate(modes): 
            c.cd()
            hist=self.hists[flag][histname]
            hist.SetLineColor(colours[i])
            if i==0: hist.Draw()
            else: hist.Draw('same')
            l.AddEntry(hist, flag)
        l.Draw()

        outfilename = '/eos/experiment/sndlhc/users/aconsnd/simulation/neutrino/data/nue-extendedreconstruction/results.root'
        outf = ROOT.TFile.Open(outfilename, 'recreate')
        outf.WriteObject(c, c.GetName(), 'kOverwrite')
        outf.Close()
        print(f'Comparison canvas written to {outfilename}')
        
def SetUpAnalysisClass():
    from argparse import ArgumentParser
    from AnalysisFunctions import Analysis 
    from args_config import add_arguments
    muAna_parser = ArgumentParser()
    add_arguments(muAna_parser)
    muAna_options = muAna_parser.parse_args()
    muAna_options.simulation=True
    muAna = Analysis(muAna_options)
    # muAna.BuildBarLengths(geo.modules['MuFilter'])
    # muAna.Makecscintdict(muAna_options.TWCorrectionRun, 'corrected')
    
    if not muAna.simulation:
        timealignment=muAna.GetTimeAlignmentType(runNr=str(muAna_options.runNr).zfill(6))
        muAna.MakeAlignmentParameterDict(timealignment)
        muAna.MakeTWCorrectionDict()
        
    return muAna       


class QuarkVectorExtrapolation(object):
    def __init__(self, options, tw):
    
        self.options=options
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

        self.hists=tw.hists

        self.sigmatds0=0.263 # ns 

        self.A, self.B = ROOT.TVector3(), ROOT.TVector3()
        self.MuFilter.GetPosition(24004, self.B, self.B)
        self.HCAL5z = 0.5*(self.A.z() + self.B.z())
        
        if options.signalpartitions: self.Loadnumuevents()

        self.numuStudy=True if options.numuStudy else False 

        # Not exact, just rough for rejecting rubbish combinations of DS clusters
        self.acceptancelimits={'x':[-100, 20], 'y':[-10, 100]}

        self.ReadInVectors()
        if not self.quark_vectors: 
            os.exit(1)

    def ReadInVectors(self):

        quarkvectorfilename = self.options.path+self.options.fname.replace('.root', '_quark_momentum.txt')
        if not os.path.exists(quarkvectorfilename):
            print(f'No struck quark vector file for {self.options.fname}')
            return 

        self.quark_vectors={}

        with open(quarkvectorfilename, 'r') as file:
            for line in file:
                # Strip any leading/trailing whitespace characters
                line = line.strip()
                
                # Split the line by commas
                event_number, key, x, y, z, px, py, pz = line.split(' ')
                
                # Convert the values to their respective types
                event_number = int(event_number)
                key = int(key)
                x = float(x)
                y = float(y)
                z = float(z)
                px = float(px)
                py = float(py)
                pz = float(pz)                
                
                # Assign to the dictionary
                self.quark_vectors[event_number] = [key, x, y, z, px, py, pz]

    def StruckQuarkExtrapolation(self, hits):

        """
        Here I will take the vector of the struck quark, evaluate the 
        vector at the z-values of the fired planes and plot the residual 
        between the barycentre and the vector. 
        """

        self.barycentres = self.muAna.GetBarycentres(hits)
        self.xbarycentres = self.muAna.GetOverallXBarycentre(self.barycentres, mode='maxQDC')

        key, pos, mom = self.GetQuarkPosMom()
        mom_mag = np.sqrt(sum([i**2 for i in mom]))
        mom_normalised = [i/mom_mag for i in mom]

        for plane in range(5):

            if any([plane not in self.barycentres, plane not in self.xbarycentres]): continue

            # Extrapolate struck quark vector to this plane
            x_q, y_q = self.ExtrapolateStruckQuark(plane, pos, mom)

            residual_x = x_q - self.xbarycentres[plane]['dxB']
            residual_y = y_q - self.barycentres[plane]['y-barycentre']['yB']

            histname = f'quarkvectorresidualx_plane{plane}'
            if not histname in self.hists:
                title = '#splitline{Residual in xz projection between struck quark momentum vector and measured barycentre}{plane '+str(plane+1)+'};Residual in xz [cm]'
                self.hists[histname] = ROOT.TH1F(histname, title, 100, -50, 50)
            self.hists[histname].Fill(residual_x)

            histname = f'quarkvectorresidualy_plane{plane}'
            if not histname in self.hists:
                title = '#splitline{Residual in yz projection between struck quark momentum vector and measured barycentre}{plane '+str(plane+1)+'};Residual in xz [cm]'
                self.hists[histname] = ROOT.TH1F(histname, title, 100, -50, 50)
            self.hists[histname].Fill(residual_y)

            histname = f'quarkvectorresidualxy_plane{plane}'
            if not histname in self.hists:
                title = '#splitline{Correlation of residuals in projections between struck quark momentum vector and measured barycentre}{plane '+str(plane+1)+'};Residual in xz [cm];Residual in yz [cm];Counts'
                self.hists[histname] = ROOT.TH2F(histname, title, 120, -60, 60, 120, -60, 60)
            self.hists[histname].Fill(residual_x, residual_y)            
            
    def ExtrapolateStruckQuark(self, plane, pos, mom):
        # Load in z pos of plane 
        self.tw.MuFilter.GetPosition(20000+1000*plane, self.A, self.B)
        z_target = 1/2 * (self.A.z() + self.B.z())
        dz = z_target - pos[2]*100

        # x_q = pos[0]*100 + scale*mom[0]
        x_q = pos[0]*100 + mom[0]/mom[2]*dz
        # y_q = pos[1]*100 + scale*mom[1] 
        y_q = pos[1]*100 + mom[1]/mom[2]*dz

        return x_q, y_q

    def GetQuarkPosMom(self):
        x = self.quark_vectors[self.tw.M.EventNumber]
        key=x[0]
        pos=x[1:4]
        mom=x[4:]
        return key, pos, mom


    def WriteOutHistograms(self):

        if not self.simulation:
            d = f'{self.outpath}splitfiles/run{self.runNr}/{self.tw.mode}/'
            os.makedirs(d, exist_ok=True)
            
            outfilename=d+f'struckquarkextr_{self.options.nStart}.root' 

        elif self.simulation and self.options.simMode=='neutrino':

            if self.options.mode=='struckquark': 

                d = f'{self.outpath}{self.tw.mode}/'
                os.makedirs(d, exist_ok=True)

                dirkey1, dirkey2, filename = self.options.fname.split('/')
                key=filename.replace('.root', '').split('_')[1]

                outfilename = d+f'struckquarkextr_{key}.root'

        elif self.simulation and self.options.simMode == 'neutralhadron':
            d = f'{self.outpath}{self.tw.mode}/'
            os.makedirs(d, exist_ok=True)

            particle_type, Emin, Emax, key = self.options.fname.split('_')[3:7]
            outfilename=d+f'struckquarkextr_{particle_type}_{Emin}_{Emax}_{key}.root'

        elif self.simulation and self.options.simMode == 'passingmuon':
            d = f'{self.outpath}{self.tw.mode}/'
            os.makedirs(d, exist_ok=True)

            keys=self.options.fname.split('_')[1:3]
            outfilename=d+f'struckquarkextr_{keys[0]}_{keys[1]}.root'

        elif self.simulation and self.options.simMode == 'nue':
            print(f'Not implemented nue saving protocol! ')     

        if os.path.exists(outfilename): outfile=ROOT.TFile.Open(outfilename, 'recreate')
        else: outfile=ROOT.TFile.Open(outfilename, 'create')            

        for hname in self.hists:
            
            hist=self.hists[hname]
            # if hname in ('yEx', 'xEx', 'reft', 'n-xz-combinations', 'n-yz-combinations', 'n_DSclusters'):
            #     outfile.WriteObject(hist, hname, 'kOverwrite')
            #     continue
            if hname.find('plane')>0:
                key, plane = hname.split('_')

                if not hasattr(outfile, key): folder=outfile.mkdir(key)
                else: folder=outfile.Get(key)
                
                folder.cd()
                hist.Write(hname, 2) # The 2 means it will overwrite a hist of the same name

        outfile.Close()
        print(f'{len(self.hists)} histograms saved to {outfilename}')    
