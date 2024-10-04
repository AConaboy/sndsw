#!/usr/bin/env python
import ROOT, csv, os, pickle
import numpy as np
from pathlib import Path

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

        self.eventswithcombinations=[]

        # Not exact, just rough for rejecting rubbish combinations of DS clusters
        self.acceptancelimits={'x':[-100, 20], 'y':[-10, 100]}

        if self.options.mode=='nue-extendedreconstruction' and self.simulation: 

            d = f'{self.outpath}{self.tw.mode}/'
            os.makedirs(d, exist_ok=True)

            dirkey1, dirkey2, filename = self.options.fname.split('/')
            key=filename.replace('.root', '').split('_')[1]

            self.keynamedict = {'wMuon':'with muon', 'woMuon': 'w/o muon', 'allEvents':'all events'}

            self.datafilename=d+f'extendedreconstruction_{key}.csv'
            column_names=['hasMuon','DSmult0x','DSmult0y','DSmult1x','DSmult1y','USmult3','USmult4','dx3','dy3','dx4','dy4','lambdax3','lambdax4', 'lambday3', 'lambday4']

            with open(self.datafilename, 'w') as f:
                writer=csv.writer(f)
                writer.writerow(column_names)

    def ExtendReconstruction(self, hits):
        # Here I want to get the points in space from the DS hits, and see if a US hit aligns with these
        # If they do then I can plot the doca between the line formed between these DS hits and the US hit
        # if self.options.OutgoingMuon=='yes' and not eventHasMuon: return
        # elif self.options.OutgoingMuon=='no' and eventHasMuon: return
        self.eventHasMuon=self.OutgoingMuon()
        if self.eventHasMuon: self.muonhistkey = 'wMuon'
        elif not self.eventHasMuon: self.muonhistkey = 'woMuon'

        self.dsClusters = self.tw.M.trackTask.clusMufi
        if len(self.dsClusters)==0: return
        
        self.barycentres = self.muAna.GetBarycentres(hits)
        self.xbarycentres = self.muAna.GetOverallXBarycentre(self.barycentres, mode='maxQDC')

        self.GetDSPoints() # Working with clusters
        fired_planes=list(self.DS_points.keys())
        if len(fired_planes)==3: 
            self.RecordEventNr()
            print(f'3 fired DS planes in event {self.tw.M.EventNumber}')
            return        
        # self.muAna.print_timestamp("Got DSpoints dict") 

        histname = 'n_DSclusters'
        if not histname in self.hists:
            title="Number of cluster formed in muon system;# DS clusters;Counts"
            self.hists[histname] = ROOT.TH1F(histname, title, 26, 0, 26)
        self.hists[histname].Fill(len(self.dsClusters))
        
        """
        For each permutation of pairs of x,y fired bars, I can make a straight line
        and extrapolate to the US to check the doca with the x-barycentre, y-barycentre.

        In order to know which permutation is successful, I will need to keep track of the detector IDs somehow. 
        """
        
        self.GetCombinatorics()

        for p in ('x', 'y'):
            histname=f'n-{p}z-combinations'
            if not histname in self.hists:
                title=f'Number of {p} clusters permutations in the DS;# {p} cluster permutations;Counts'
                self.hists[histname] = ROOT.TH1F(histname, title, 6, 0, 6)
            self.hists[histname].Fill(len(self.combinations[p]))

        if len(self.combinations)>0: self.RecordEventNr()
        else: return 

        # Count HCAL hits in each plane
        self.GetMultiplicity(hits)

        # If 2 fired planes in the DS, I can connect the space points in each combination together and look for a hit in the US
        
        if len(fired_planes)==2:
            self.xy_residuals = {plane:{'x':False, 'y':False} for plane in range(5)}

            for idx, proj in enumerate(['x', 'y']):
                
                residuals=[]
                for combination in self.combinations[proj]:
                    
                    # Make lines for the xz and yz projections that join the points of this pair
                    line = self.ConnectPoints(combination, proj)
                    # line == False if the points draw a line out of the acceptance of HCAL plane 5
                    if not line: continue

                    # returns a dictionary of the residual in that projection
                    res = self.USresidual(line, proj) 
                    if not 4 in res and 5 in res: continue

                    residuals.append(res)

                if len(residuals)==0: continue

                # Define best combination in each projection as the one with the lowest residual. Projections are orthogonal so no issue.
                best_residual  = min(residuals, key=lambda d:d[4])

                for plane in best_residual:
                    # Update value for each projection if there are suitable combinations
                    self.xy_residuals[plane][proj] = best_residual[plane]

            self.lambda_x_dict = {i:None for i in range(5)}
            self.lambda_y_dict = {i:None for i in range(5)}

            for plane in self.xy_residuals:

                lambda_x = abs(self.xbarycentres[plane]['lambda_x'])
                self.lambda_x_dict[plane] = lambda_x
                lambda_y = self.barycentres[plane]['y-barycentre']['lambda_y']
                self.lambda_y_dict[plane] = lambda_y

                if len(self.xy_residuals[plane])==2:

                    self.xyresiduals_hists(plane)
                    # self.USmultvrr_hists(plane)

                self.lambda_hists(plane)

                self.ResultantResidual_hists(plane)

            self.writedata()

        elif len(fired_planes)==1:
            pass

    def writedata(self):

        output_line = self.filekey,self.tw.M.EventNumber,self.eventHasMuon,self.multiplicity_dict[3][0],self.multiplicity_dict[3][1],self.multiplicity_dict[3][2],self.multiplicity_dict[3][3],self.multiplicity_dict[2][3],self.multiplicity_dict[2][4],self.xy_residuals[3]['x'],self.xy_residuals[3]['y'],self.xy_residuals[4]['x'],self.xy_residuals[4]['y'],self.lambda_x_dict[3], self.lambda_x_dict[4],self.lambda_y_dict[3], self.lambda_y_dict[4]

        with open(self.datafilename, 'a', newline='') as f:
            writer=csv.writer(f)
            writer.writerow(output_line)

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
        self.hists[histname].Fill(*self.xy_residuals[plane].values())                
        self.hists[f'xyresidual_allEvents_plane{plane}'].Fill(*self.xy_residuals[plane].values())

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
        lambda_x = self.lambda_x_dict[plane]
        if lambda_x: 
            self.hists[histname].Fill(lambda_x)
            self.hists[f'lambdax_allEvents_plane{plane}'].Fill(lambda_x)
        
        histname = f'lambday_{self.muonhistkey}_plane{plane}'
        lambda_y = self.lambda_y_dict[plane]
        if lambda_y: 
            self.hists[histname].Fill(lambda_y)
            self.hists[f'lambday_allEvents_plane{plane}'].Fill(lambda_y)

    def ResultantResidual_hists(self, plane):

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
        
        if 'x' in self.xy_residuals[plane]:
            histname = f'dx_{self.muonhistkey}_plane{plane}'
            self.hists[histname].Fill(self.xy_residuals[plane]['x'])
            self.hists[f'dx_allEvents_plane{plane}'].Fill(self.xy_residuals[plane]['x'])
        
        if 'y' in self.xy_residuals[plane]:
            histname = f'dy_{self.muonhistkey}_plane{plane}'
            self.hists[histname].Fill(self.xy_residuals[plane]['y'])
            self.hists[f'dy_allEvents_plane{plane}'].Fill(self.xy_residuals[plane]['y'])
        
        if 'x' in self.xy_residuals[plane] and 'y' in self.xy_residuals[plane]:
            dx, dy = self.xy_residuals[plane]['x'], self.xy_residuals[plane]['y']
            ds = np.sqrt(dx**2 + dy**2)
            histname = f'ds_{self.muonhistkey}_plane{plane}'
            self.hists[histname].Fill(ds)
            self.hists[f'ds_allEvents_plane{plane}'].Fill(ds)                          

    def InAcceptance(self, line, proj):
        
        if any({
            line(self.HCAL5z) < self.acceptancelimits[proj][0],
            line(self.HCAL5z) > self.acceptancelimits[proj][1],
                }): return False 
        else: return True 

    def hasDShits(self, hits):
        
        detIDs = [hit.GetDetectorID() for hit in hits if hit.isValid()]
        nDShits = 0
        
        for i in detIDs: 
            if self.muAna.parseDetID(i)[0]==3: nDShits+=1
    
        return nDShits

    def GetDSPoints(self):
        
        self.DS_points = {}

        for c in self.dsClusters:
            first=c.GetFirst()
            s,p,b = self.muAna.parseDetID(first)

            if not p in self.DS_points: self.DS_points[p] = {'x':[], 'y':[]}

            # self.MuFilter.GetPosition(first, self.A, self.B)
            c.GetPosition(self.A, self.B)
            avg_z = 0.5*(self.A.z() + self.B.z())

            if b<60: 
                avg_y = 0.5*(self.A.y() + self.B.y())
                self.DS_points[p]['y'].append( [avg_y, avg_z] )
                # self.DS_points[p]['y'].append(avg_y)
            elif b>=60: 
                avg_x = 0.5*(self.A.x() + self.B.x())
                self.DS_points[p]['x'].append( [avg_x, avg_z] )            

    def NDShits(self, hits):
        dshits = [i for i in hits if self.muAna.parseDetID(i.GetDetectorID())[0]==3]

    def GetCombinatorics(self):
        
        fired_planes=list(self.DS_points.keys())

        xz_proj, yz_proj = self.DS_points[fired_planes[0]]['x'], self.DS_points[fired_planes[0]]['y']

        # Combinatorics are just different within 1 plane
        if len(fired_planes)==1: 
            
            """
            Here all I can do is return the cluster centres in the fired DS plane
            """
            self.combinations = {'x': xz_proj, 'y': yz_proj}

        elif len(fired_planes)==2:
            """
            For 2 fired DS planes, each horizontal bar can pair with each vertical bar. 
            For each fired horizontal bar, there are as many pairs as there are fired vertical bars in the same station. 

            For each pair of fired horizontal and vertical bars in plane i, there are N_nextplane combinations with a pair in plane i+1
            Where N_nextplane is the number of horizontal and vertical bar pairs that can be formed in the next plane

            """

            xz_proj_next, yz_proj_next = self.DS_points[fired_planes[1]]['x'], self.DS_points[fired_planes[1]]['y']

            xz_combs = [[xz_val, xz_val_next] for xz_val in xz_proj for xz_val_next in xz_proj_next]
            yz_combs = [[yz_val, yz_val_next] for yz_val in yz_proj for yz_val_next in yz_proj_next]

            self.combinations = {'x': xz_combs, 'y': yz_combs}
            
    def ConnectPoints(self, combination, proj):
        
        """
        Each combination passed to this method is a pair of points 
        for successive planes in the DS. 

        The structure of combination here is: 
        A,B,C,D = [ [[xi,zi], [yi,zi]], [[xj,zj], [yj,zj]] ] 

        So the xz proj line is between: A,C 
        And the yz proj line is between: B,D
        """
        
        plane0, plane1 = combination
        
        # Make line for both projections

        m = (plane1[0] - plane0[0])/ (plane1[1] - plane0[1]) 
        c = plane1[0] - (m * plane1[1])
        line = lambda z : m*z+c

        if not self.InAcceptance(line, proj): return False
        else: return line

    def USresidual(self, line, proj):

        res_dict={}
        for plane in range(5):
    
            # Get barycentre 
            b=self.GetBarycentre(plane, proj)
            if not b:continue

            # Get z position of plane to pass into eqn of line
            self.tw.MuFilter.GetPosition(20000+1000*plane, self.A, self.B)
            HCAL_z = 0.5*(self.A.z() + self.B.z())
            
            ext = line(HCAL_z) 
            residual = ext-b
            res_dict[plane]=residual
            
        return res_dict

    def GetBarycentre(self, plane, proj):

        if proj=='x': 
            if plane not in self.xbarycentres: return
            if not 'dxB' in self.xbarycentres[plane]:return
            b=self.xbarycentres[plane]['dxB']
        elif proj=='y': 
            if plane not in self.barycentres: return
            if not 'y-barycentre' in self.barycentres[plane]:return
            if not 'yB' in self.barycentres[plane]['y-barycentre']:return
            b=self.barycentres[plane]['y-barycentre']['yB']
        return b

    def GetLambda(self, plane, proj):
        if proj=='x': 
            if plane not in self.xbarycentres: return
            if not 'lambda_x' in self.xbarycentres[plane]:return
            b=self.xbarycentres[plane]['lambda_x']
        elif proj=='y': 
            if plane not in self.barycentres: return
            if not 'y-barycentre' in self.barycentres[plane]:return
            if not 'lambda_y' in self.barycentres[plane]['y-barycentre']:return
            b=self.barycentres[plane]['y-barycentre']['lambda_y']
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
                
        if self.eventHasMuon: self.hists[f'USmultiplicity_wMuon_plane{plane}'].Fill(self.multiplicity_dict[2][plane])
        elif not self.eventHasMuon: self.hists[f'USmultiplicity_woMuon_plane{plane}'].Fill(self.multiplicity_dict[2][plane])
        self.hists[f'USmultiplicity_allEvents_plane{plane}'].Fill(self.multiplicity_dict[2][plane])
        
    def RecordEventNr(self):
        fired_planes=list(self.DS_points.keys())
        event_data = [self.options.fname, self.tw.M.EventNumber, len(fired_planes)]
        self.eventswithcombinations.append(event_data)

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

        print(f'Data written to {self.datafilename}') 

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

        self.eventswithcombinations=[]

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