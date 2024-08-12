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
        
        if options.signalpartitions: self.Loadnumuevents()
        
        self.numuStudy=True if options.numuStudy else False 

    def ExtendReconstruction(self, hits):
        # Here I want to get the points in space from the DS hits, and see if a US hit aligns with these
        # If they do then I can plot the doca between the line formed between these DS hits and the US hit

        # self.muAna.print_timestamp("Starting emr method")    

        self.dsClusters = self.tw.M.trackTask.clusMufi
        if len(self.dsClusters)==0: return

        # nDShits = self.hasDShits(hits)
        # if nDShits==0: return 
        
        self.barycentres = self.muAna.GetBarycentres(hits)
        # self.muAna.print_timestamp("Got barycentres") 

        self.GetDSPoints() # Working with clusters
        # self.muAna.print_timestamp("Got DSpoints dict") 
        """
        For each permutation of pairs of x,y fired bars, I can make a straight line
        and extrapolate to the US to check the doca with the x-barycentre, y-barycentre.

        In order to know which permutation is successful, I will need to keep track of the detector IDs somehow. 
        """
        self.GetCombinatorics()
        print(f"len(combinatorics): {len(self.combinations)}")
        # self.muAna.print_timestamp("Got combinatorics") 

        # If 2 fired planes in the DS, I can connect the space points in each combination together and look for a hit in the US
        fired_planes=list(self.DS_points.keys())
        if fired_planes==2:
            for combination in self.combinations:
                # Make lines for the xz and yz projections that join the points of this pair
                lines = self.ConnectPoints(combination)
                doca = self.USdoca(hits, lines)


    def hasDShits(self, hits):
        
        detIDs = [hit.GetDetectorID() for hit in hits if hit.isValid()]
        nDShits = 0
        # [nDShits+=1 for i in detIDs if self.muAna.parseDetID(detID)[0]==3]
        for i in detIDs: 
            if self.muAna.parseDetID(i)[0]==3: nDShits+=1
    
        return nDShits

    def GetDSPoints(self):
        
        self.DS_points = {}

        for c in self.dsClusters:
            first=c.GetFirst()
            s,p,b = self.muAna.parseDetID(first)

            if not p in self.DS_points: self.DS_points = { p:{'x':[], 'y':[]} }            

            self.MuFilter.GetPosition(first, self.A, self.B)
            avg_z = 0.5*(self.A.z() + self.B.z())
            if p<60: 
                avg_y = 0.5*(self.A.y() + self.B.y())
                self.DS_points[p]['y'].append( [avg_y, avg_z] )
                # self.DS_points[p]['y'].append(avg_y)
            elif p>=60: 
                avg_x = 0.5*(self.A.x() + self.B.x())
                self.DS_points[p]['x'].append( [avg_x, avg_z] )            
        
        """
        ### The code below is for working with individual hits

        for hit in hits:

            if not hit.isValid(): continue
            detID = hit.GetDetectorID()    
            s,p,b=self.muAna.parseDetID(detID)
            if not s==3: continue 
            
            # Storing data for plane p 
            if not p in self.DS_points: self.DS_points = { p:{'x':[], 'y':[]} }
            
            self.MuFilter.GetPosition(detID, self.A, self.B)
            avg_z = 0.5*(self.A.z() + self.B.z())
            if p<60: 
                avg_y = 0.5*(self.A.y() + self.B.y())
                self.DS_points[p]['y'].append( [avg_y, avg_z] )
                # self.DS_points[p]['y'].append(avg_y)
            elif p>=60: 
                avg_x = 0.5*(self.A.x() + self.B.x())
                self.DS_points[p]['x'].append( [avg_x, avg_z] )
                # self.DS_points[p]['x'].append(avg_x)
        """
        """
        If a projection in a DS station has multiple bars firing, 
        I should store sufficient info on these bars in order to test the combinatorics
        """
        
    def GetCombinatorics(self):
        
        fired_planes=list(self.DS_points.keys())

        combinations=[]

        # x_proj, y_proj = [i[0] for i in self.DS_points[fired_planes[0]]['x']], [i[0] for i in self.DS_points[fired_planes[0]]['y']]
        xz_proj, yz_proj = self.DS_points[fired_planes[0]]['x'], self.DS_points[fired_planes[0]]['y']

        # Combinatorics are just different within 1 plane
        if fired_planes==1: 
            [combinations.append([xz_val, yz_val]) for xz_val in xz_proj for yz_val in yz_proj]

        elif fired_planes==2:

            """
            For 2 fired DS planes, each horizontal bar can pair with each vertical bar. 
            For each fired horizontal bar, there are as many pairs as there are fired vertical bars in the same station. 

            For each pair of fired horizontal and vertical bars in plane i, there are N_nextplane combinations with a pair in plane i+1
            Where N_nextplane is the number of horizontal and vertical bar pairs that can be formed in the next plane

            """

            xz_proj_next, yz_proj_next = self.DS_points[fired_planes[1]]['x'], self.DS_points[fired_planes[1]]['y']

            next_plane_combs = [[xz_val_next, yz_val_next] for xz_val_next in xz_proj_next for yz_val_next in yz_proj_next]
            plane_combs = [[xz_val, yz_val] for xz_val in xz_proj for yz_val in yz_proj]

            [combinations.append([i, j]) for i in plane_combs for j in next_plane_combs]

        self.combinations = combinations

    def ConnectPoints(self, combination):
        
        """
        Each combination passed to this method is a pair of points 
        for successive planes in the DS. 

        The structure of combination here is: 
        A,B,C,D = [ [[xi,zi], [yi,zi]], [[xj,zj], [yj,zj]] ] 

        So the xz proj line is between: A,C 
        And the yz proj line is between: B,D
        """
        
        plane0, plane1 = combination
        
        plane0_xz, plane0_yz = plane0
        plane1_xz, plane1_yz = plane1

        # Make line for both projections

        m_xz = (plane1_xz[0] - plane0_xz[0])/ (plane1_xz[0] - plane0_xz[0]) 
        c_xz = plane1_xz[0] / (m_xz * plane0_xz[1])
        line_xz = lambda x : m_xz*x+c_xz 

        m_yz = (plane1_yz[0] - plane0_yz[0])/ (plane1_yz[0] - plane0_yz[0])
        c_yz = plane1_yz[0] / (m_yz * plane0_yz[1])
        line_yz = lambda x : m_yz*x+c_yz 

        return {'xz':line_xz, 'yz':line_yz}

    def USdoca(self, hits, lines): 

        docas = {} 

        for plane in range(5):
            xB = self.barycentres[plane][]
        