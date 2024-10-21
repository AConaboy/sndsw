#!/usr/bin/env python
import ROOT, os, csv, json, joblib
from datetime import datetime
import math as m 
import numpy as np
from itertools import combinations
import pandas as pd


ROOT.gInterpreter.ProcessLine('#include "/afs/cern.ch/user/a/aconsnd/sndsw/analysis/tools/sndSciFiTools.h"')

class HCALTools(object):
    def __init__(self,muAna,MuFilter):
        self.muAna = muAna
        self.A, self.B = ROOT.TVector3(), ROOT.TVector3()
        self.clusMufi=ROOT.TObjArray(100)
        self.simulation=muAna.simulation

        MuFilter.GetPosition(24004, self.B, self.B)
        self.HCAL5z = 0.5*(self.A.z() + self.B.z())
        self.acceptancelimits={'x':[-100, 20], 'y':[-10, 100]}

        self.column_names=['filekey', 'EventNumber','interactionWall',
        'scifi_median_x','scifi_median_y',
        'scifi_residual_x','scifi_residual_y',
        'dx0','dx1','dx2','dx3','dx4',
        'dy0','dy1','dy2','dy3','dy4',
        'x0','x1','x2','x3','x4',
        'y0','y1','y2','y3','y4',
        'lambdax0','lambdax1', 'lambdax2','lambdax3', 'lambdax4', 
        'lambday0','lambday1', 'lambday2','lambday3', 'lambday4',
        'HCAL5barcode'
        ]

        self.muAna.BuildBarLengths(MuFilter)
        self.muAna.Makecscintdict('005408')

    def setsimulation(self, sim, run=-1):
        self.muAna.simulation=sim
        self.simulation=sim

        if sim==False:
            if run==-1:
                print(f'Have to pass a valid runNr for real data!')
                return
            if not hasattr(self.muAna, 'twparameters'): self.muAna.MakeTWCorrectionDict()
            if not hasattr(self.muAna, 'alignmentparameters'): 
                self.muAna.timealignment=self.muAna.GetTimeAlignmentType(runNr=str(run).zfill(6))
                self.muAna.MakeAlignmentParameterDict()

    def GetDSClusterCentroids(self):
        
        self.DS_centroids = {}

        for c in self.clusMufi:
            first=c.GetFirst()
            s,p,b = self.muAna.parseDetID(first)

            if not p in self.DS_centroids: self.DS_centroids[p] = {'x':[], 'y':[]}

            # self.MuFilter.GetPosition(first, self.A, self.B)
            c.GetPosition(self.A, self.B)
            avg_z = 0.5*(self.A.z() + self.B.z())

            if b<60: 
                avg_y = 0.5*(self.A.y() + self.B.y())
                self.DS_centroids[p]['y'].append( [avg_y, avg_z] )
            elif b>=60: 
                avg_x = 0.5*(self.A.x() + self.B.x())
                self.DS_centroids[p]['x'].append( [avg_x, avg_z] )

    def GetCombinatorics(self):
        
        fired_planes=list(self.DS_centroids.keys())

        xz_proj, yz_proj = self.DS_centroids[fired_planes[0]]['x'], self.DS_centroids[fired_planes[0]]['y']

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

            xz_proj_next, yz_proj_next = self.DS_centroids[fired_planes[1]]['x'], self.DS_centroids[fired_planes[1]]['y']

            xz_combs = [[xz_val, xz_val_next] for xz_val in xz_proj for xz_val_next in xz_proj_next]
            yz_combs = [[yz_val, yz_val_next] for yz_val in yz_proj for yz_val_next in yz_proj_next]

            self.combinations = {'x': xz_combs, 'y': yz_combs}   

    def Get_interaction_median_positions(self, scifi_hits, Scifi):
        # scifi=self.tw.Scifi

        if self.interactionWall==6: wall=5
        else: wall=self.interactionWall
        
        scifi_hitDict = {i.GetDetectorID():i for i in scifi_hits if i.GetDetectorID()//1000000==wall} 

        positions = {'x':[], 'y':[]}
        for detID, hit in scifi_hitDict.items():
            Scifi.GetSiPMPosition(detID, self.A,self.B)
            if hit.isVertical():
                pos = 0.5*(self.A[0] + self.B[0])
                positions['x'].append(pos)
            else:
                pos = 0.5*(self.A[1] + self.B[1]) 
                positions['y'].append(pos)
        
        self.intWall_median_x, self.intWall_median_y = np.median(positions['x']),np.median(positions['y'])  


    def BDT_cut(self, event, MuFilter, Scifi):

        hits = event.Digi_MuFilterHits
        scifi_hits = event.Digi_ScifiHits

        # Get event t0 time if not simulation
        if not self.simulation:
            if self.muAna.options.referencesystem==1:
                self.reft = self.muAna.GetScifiAverageTime(Scifi, scifi_hits)
            elif self.muAna.options.referencesystem==3:
                self.reft = self.muAna.GetDSHaverage(hits)

        # Get predicted interaction wall
        # Returns number between 0,4 or 5 for unreconstructed interaction wall

        filtered_hits = ROOT.snd.analysis_tools.filterScifiHits(scifi_hits, 0, "TI18")
        self.interactionWall = ROOT.snd.analysis_tools.showerInteractionWall(filtered_hits, 0, "TI18") + 1 
        
        # Get the median scifi position in x and y for the interaction wall
        self.Get_interaction_median_positions(scifi_hits=scifi_hits, Scifi=Scifi)

        # Get the lambda values and barycentres for all HCAL planes
        self.barycentres=self.muAna.GetBarycentres(hits, MuFilter=MuFilter)
        self.xbarycentres=self.muAna.GetOverallXBarycentre(self.barycentres, mode='maxQDC')

        # Get DS cluster centroids
        self.dsCluster(hits, MuFilter)
        self.GetDSClusterCentroids()
        fired_planes=list(self.DS_centroids.keys())

        # Can't use BDT here, must assume True
        if len(fired_planes)==0:
            return True   

        # Sanity check, should be none
        if len(fired_planes)==3:
            print(f'3 fired DS planes in event {self.EventNumber}')
            return True

        # Find combinations of DS cluster centroids
        self.GetCombinatorics()

        if len(fired_planes)==2:
            self.xy_residuals = {plane:{'x':np.nan, 'y':np.nan} for plane in range(5)}
            self.xy_residuals['intWall']={'x':np.nan, 'y':np.nan}

            for idx, proj in enumerate(['x', 'y']):
                
                residuals=[]
                for combination in self.combinations[proj]:
                    
                    # Make lines for the xz and yz projections that join the points of this pair
                    line = self.ConnectPoints(combination, proj)
                    # line == False if the points draw a line out of the acceptance of HCAL plane 5
                    if not line: continue

                    # returns a dictionary of the residual in that projection
                    res = self.USresidual(line, proj, MuFilter, Scifi) 
                    if not 4 in res and 5 in res: continue

                    residuals.append(res)

                if len(residuals)==0: continue

                # Define best combination in each projection as the one with the lowest residual. Projections are orthogonal so no issue.
                best_residual  = min(residuals, key=lambda d:d['intWall'])

                for plane in best_residual:
                    # if plane=='intWall':continue
                    # Update value for each projection if there are suitable combinations
                    self.xy_residuals[plane][proj] = best_residual[plane]
 
            # Require that the xy_residual is defined for the 4th and 5th plane
            if list(self.xy_residuals[4].values()) == [np.nan, np.nan]: return True
            if list(self.xy_residuals[3].values()) == [np.nan, np.nan]: return True 

            self.lambda_x_dict = {i:np.nan for i in range(5)}
            self.lambda_y_dict = {i:np.nan for i in range(5)}

            for plane in self.xy_residuals:

                if plane=='intWall':continue
                if not plane in self.xbarycentres: continue
                if not plane in self.barycentres: continue
                
                lambda_x = abs(self.xbarycentres[plane]['lambda_x'])
                self.lambda_x_dict[plane] = lambda_x
                lambda_y = self.barycentres[plane]['y-barycentre']['lambda_y']
                self.lambda_y_dict[plane] = lambda_y

                # if len(self.xy_residuals[plane])==2:

                    # self.xyresiduals_hists(plane)
                    # self.USmultvds_hists(plane)

                # if mode=='write': 
                #     self.lambda_hists(plane)
                #     self.ds_hists(plane)

            self.GetHCAL5barscode(hits)
            self.Get_ds()

            x=self.getdata(mode='get')
            xdf = pd.DataFrame([x])
            xdf['HCAL5barcode'] = xdf['HCAL5barcode'].astype(int)
            
            cols2drop = ['filekey', 'EventNumber']
            if self.simulation: cols2drop.append('hasMuon')
            xdf.drop(cols2drop, axis=1, inplace=True)
            
            features = self.model.feature_names_in_
            xdf = xdf[features] # Ensure data is in the right order for BDT 

            res = self.model.predict(xdf)
            return bool(res[0])
            
        # Also must assume True here at the moment.
        # Will test extending BDT to 1 fired DS plane
        elif len(fired_planes)==1:
            return True
        else: pass        

    def Get_ds(self):
        self.ds = dict.fromkeys([f'ds{i}' for i in range(5)], np.nan)
        for plane in range(5):
            self.ds[f'ds{plane}'] = np.sqrt( sum([i**2 for i in self.xy_residuals[plane].values()]) )

        self.scifi_residual = np.sqrt( sum([i**2 for i in self.xy_residuals['intWall'].values()]) )

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

    def USresidual(self, line, proj, MuFilter, Scifi):

        res_dict={}

        for plane in range(5):
    
            # Get barycentre 
            b=self.GetBarycentre(plane, proj)
            if not b:continue

            # Get z position of plane to pass into eqn of line
            MuFilter.GetPosition(20000+1000*plane, self.A, self.B)
            HCAL_z = 0.5*(self.A.z() + self.B.z())
            
            ext = line(HCAL_z) 
            residual = ext-b
            res_dict[plane]=residual
        
        if self.interactionWall==6: wall=5
        else: wall=self.interactionWall

        dummy_detID=int(1e6*wall + 1e5*0 + 1e4 + 1e3 + 1)
        print(f'Getting Scifi position for detID {dummy_detID}')
        Scifi.GetPosition(dummy_detID, self.A, self.B)
        scifi_intwall_z = 0.5*(self.A.z() + self.B.z())
        
        ext = line(scifi_intwall_z)
        if proj=='x': scifi_b = self.intWall_median_x
        elif proj=='y': scifi_b = self.intWall_median_y
        
        residual = ext-scifi_b
        res_dict['intWall']=residual

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

    def GetHCAL5barscode(self, hits):
        c=[False]*10 
        US_detID_list = [i.GetDetectorID() for i in hits if all([self.muAna.parseDetID(i.GetDetectorID())[0]==2, self.muAna.parseDetID(i.GetDetectorID())[1]==4])]
        for detID in US_detID_list:
            bar=self.muAna.parseDetID(detID)[2]
            c[bar]=True 
        
        self.HCAL5barscode=''.join(['1' if state else '0' for state in c])

    def InAcceptance(self, line, proj):
        
        if any({
            line(self.HCAL5z) < self.acceptancelimits[proj][0],
            line(self.HCAL5z) > self.acceptancelimits[proj][1],
                }): return False 
        else: return True 

    def getdata(self, mode='write'):

        x_list = [round(self.xbarycentres[p]['dxB'],3) for p in range(5)]
        y_list = [round(self.barycentres[p]['y-barycentre']['yB'],3) for p in range(5)]
        dx_list = [round(self.xy_residuals[p]['x'], 3) for p in range(5)]
        dy_list = [round(self.xy_residuals[p]['y'], 3) for p in range(5)]
        ds_list = [round(self.ds[f'ds{p}'], 3) for p in range(5)]

        lambda_x_list = [round(self.lambda_x_dict[p], 3) for p in range(5)]
        lambda_y_list = [round(self.lambda_y_dict[p], 3) for p in range(5)]

        if self.simulation:
            
            output_dict = dict.fromkeys(self.column_names, None)
            output_dict['filekey'] = self.filekey
            output_dict['EventNumber'] = self.EventNumber
            output_dict['hasMuon'] = self.eventHasMuon
            output_dict['interactionWall'] = self.interactionWall
            output_dict['scifi_median_x'] = round(self.intWall_median_x,3)
            output_dict['scifi_median_y'] = round(self.intWall_median_y,3)
            
            output_dict['ds_scifi'] = round(self.scifi_residual, 3)
            output_dict['scifi_residual_x'] = round(self.xy_residuals['intWall']['x'], 3)
            output_dict['scifi_residual_y'] = round(self.xy_residuals['intWall']['y'], 3)

            for i in range(5):
                output_dict[f'dx{i}'] = dx_list[i]
                output_dict[f'dy{i}'] = dy_list[i]
                output_dict[f'ds{i}'] = ds_list[i]
                output_dict[f'x{i}'] = x_list[i] 
                output_dict[f'y{i}'] = y_list[i] 
                output_dict[f'lambdax{i}'] = round(lambda_x_list[i], 3)
                output_dict[f'lambday{i}'] = round(lambda_y_list[i], 3)

            output_dict['HCAL5barcode'] = self.HCAL5barscode
            
            # output_line = self.filekey,self.tw.M.EventNumber,self.eventHasMuon,self.interactionWall,round(self.intWall_median_x,3), round(self.intWall_median_y,3),*scifi_residual,*x_residuals_list,*y_residuals_list,*xbarycentre_list,*ybarycentre_list,*[round(i,3) for i in self.lambda_x_dict.values()],*[round(i,3) for i in self.lambda_y_dict.values()],self.HCAL5barscode

        else:

            output_dict = dict.fromkeys(self.column_names, None)
            output_dict['filekey'] = self.filekey
            output_dict['EventNumber'] = self.EventNumber
            # output_dict['hasMuon'] = self.eventHasMuon
            output_dict['interactionWall'] = self.interactionWall
            output_dict['scifi_median_x'] = round(self.intWall_median_x,3)
            output_dict['scifi_median_y'] = round(self.intWall_median_y,3)
            
            output_dict['ds_scifi'] = round(self.scifi_residual, 3)
            output_dict['scifi_residual_x'] = round(self.xy_residuals['intWall']['x'], 3)
            output_dict['scifi_residual_y'] = round(self.xy_residuals['intWall']['y'], 3)

            for i in range(5):
                output_dict[f'dx{i}'] = dx_list[i]
                output_dict[f'dy{i}'] = dy_list[i]
                output_dict[f'ds{i}'] = ds_list[i]
                output_dict[f'x{i}'] = x_list[i] 
                output_dict[f'y{i}'] = y_list[i] 
                output_dict[f'lambdax{i}'] = round(lambda_x_list[i], 3)
                output_dict[f'lambday{i}'] = round(lambda_y_list[i], 3)

            output_dict['HCAL5barcode'] = self.HCAL5barscode

            # output_line = self.filekey,self.tw.M.EventNumber,self.interactionWall,round(self.intWall_median_x,3), round(self.intWall_median_y,3),*scifi_residual,*x_residuals_list,*y_residuals_list,*xbarycentre_list,*ybarycentre_list,*[round(i,3) for i in self.lambda_x_dict.values()],*[round(i,3) for i in self.lambda_y_dict.values()],self.HCAL5barscode

        if mode=='write':
            with open(self.datafilename, 'a', newline='') as f:
                writer=csv.writer(f)
                writer.writerow(output_dict.values())
        elif mode=='get':
            return output_dict

    def OutgoingMuon(self, event):
        ### Return true for events with an outgoing muon
        if event.GetName()=='rawConv':
            print(f'No MCTrack branch in real data! ')
            return False
        if abs(event.MCTrack[1].GetPdgCode()) == 13: return True 
        else: return False            

    def dsCluster(self, hits, mufi):
        clusters = []
        hitDict = {}
        for k in range(hits.GetEntries()):
            d = hits[k]
            # if (d.GetDetectorID()//10000)<3 or (not d.isValid()): continue
            if (d.GetDetectorID()//10000)<3: continue
            hitDict[d.GetDetectorID()] = k
        hitList = list(hitDict.keys())
        if len(hitList)>0:
            hitList.sort() # Smallest detID first: sorted low z to high z.
            tmp = [ hitList[0] ] 
            cprev = hitList[0]
            ncl = 0
            last = len(hitList)-1
            hitvector = ROOT.std.vector("MuFilterHit*")()
            for i in range(len(hitList)):
                if i==0 and len(hitList)>1: continue # tmp is the 1st hit, so start from the 2nd
                c=hitList[i] # position of hit in the Digi_MuFilterHits array (k in 1st loop)
                neighbour = False
                if (c-cprev)==1 or (c-cprev)==2:    
                    neighbour = True    # allow for one missing channel to be classed as neighbour
                    tmp.append(c) 
                if not neighbour  or c==hitList[last] or c%1000==59: # if not neighbour, if last horizontal scintillator in plane
                    first = tmp[0] 
                    N = len(tmp)
                    hitvector.clear()
                    for aHit in tmp: hitvector.push_back( hits[hitDict[aHit]])
                    aCluster = ROOT.sndCluster(first,N,hitvector,mufi,False)
                    clusters.append(aCluster)
                    if c!=hitList[last]:
                        ncl+=1
                        tmp = [c]
                    elif not neighbour:   # save last channel
                        hitvector.clear()
                        hitvector.push_back( hits[hitDict[c]])
                        aCluster = ROOT.sndCluster(c,1,hitvector,mufi,False)
                        clusters.append(aCluster)
                cprev = c
        self.clusMufi.Delete()
        for c in clusters:  self.clusMufi.Add(c)