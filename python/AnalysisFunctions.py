import ROOT,os,csv, math

# class mufi_analysis: # I could make several classes for generic functions like parseDetID and more specific ones for 

class Analysis(object):

   def __init__(self, options):
      self.options=options
      self.runNr = str(options.runNumber).zfill(6)

      afswork='/afs/cern.ch/work/a/aconsnd/Timing'
      afsuser='/afs/cern.ch/user/a/aconsnd/twfiles'
      # if self.options.path.find('commissioning/TI18')>0:
      if options.datalocation=='commissioning':
         self.path=afswork+'-commissioning/'
         # self.path='TI18'
      elif options.datalocation=='physics': 
         self.path=afswork+'-physics2022/'
         # self.path='TI18'
      elif options.datalocation=='H8':
         self.path=afswork+'-H8/'
         # self.path='H8'      

   	# mapping SiPM channels to the whole system

      self.A, self.B = ROOT.TVector3(), ROOT.TVector3()
      self.systemAndPlanes = {1:2,2:5,3:7}
      self.systemAndBars = {1:7,2:10,3:60}
      self.systemAndChannels = {1:[0,8],2:[2,6],3:[0,1]}
      self.systemAndSiPMs={1:range(16),2:(0,1,3,4,6,7,8,9,11,12,14,15),3:(1,)}
      self.verticalBarDict={0:1, 1:3, 2:5, 3:6}
      self.gelsides={0:'R', 1:'L', 2:'R', 3:'L', 4:'L'}

      self.verbose=False

   def DSHcheck(self, detID): # True if detID is a horizontal bar
      s,p,b=self.parseDetID(detID)
      if s!=3: return False 
      if p<3 and b<60:  return True
      # elif p==3 and b<60: return True
      else: return False

   def DSVcheck(self, detID): # True if detID is a vertical bar
      s,p,b=self.parseDetID(detID)
      if s!=3: return False
      if p<3 and b>59: return True 
      elif p==3: return True 
      else: return False
   
   def GetListOfChannels(self, subsystem): #Only returns horizontal channels for the DS (s==3)
      channels=[f'{self.MakeFixedCh((subsystem, plane, bar, SiPM))}' for plane in range(self.systemAndPlanes[subsystem]) for bar in range(self.systemAndBars[subsystem]) for SiPM in self.systemAndSiPMs[subsystem] ]
      return channels
   
   def BuildBarLengths(MuFilter):
      Vetobarlength = MuFilter.GetConfParF('MuFilter/VetoBarX')
      USbarlength = MuFilter.GetConfParF('MuFilter/UpstreamBarX')
      DSbarlength_hor = MuFilter.GetConfParF('MuFilter/DownstreamBarX')
      DSbarlength_vert = MuFilter.GetConfParF('MuFilter/UpstreamBarY_ver')
      barlengths={1:Vetobarlength, 2:USbarlength, 3:DSbarlength_hor, 4:DSbarlength_vert}
      return barlengths

   def BuildzPos(self, MuFilter, Scifi):
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
             MuFilter.GetPosition(s*10000+p*1000+bar,A,B)
             zPos['MuFilter'][s*10+plane] = (self.A.Z()+self.B.Z())/2.
      for s in range(1,6):
         mat   = 2
         sipm = 1
         channel = 64
         for o in range(2):
             Scifi.GetPosition(channel+1000*sipm+10000*mat+100000*o+1000000*s,A,B)
             zPos['Scifi'][s*10+o] = (A.Z()+B.Z())/2.
      return zPos

   def IsSmallSiPMchannel(self, i):
      if i==2 or i==5 or i==10 or i==13: return True
      else: return False

   def GetSide(self, fixed_ch):
      detID, SiPM = fixed_ch.split('_')
      s,p,b=self.parseDetID(int(detID))
      if s==2:
         if int(SiPM)<8: side='left'
         elif int(SiPM)>7: side='right'
      elif s==1:
         if int(SiPM)<6: side='left'
         elif int(SiPM)>5: side='right'
      elif s==3:
         s,p,b=self.parseDetID(int(detID))
         if p!=3 and b<60 and int(SiPM)==0: side='left'
         elif p!=3 and b<60 and int(SiPM)==1: side='right'
         elif p!=3 and b>59 and int(SiPM)==0: side='top'
         elif p==3 and int(SiPM)==0: side='top'
         else: print('huh?')

      return side

   def IsGel(self, fixed_ch):
      detID, SiPM = fixed_ch.split('_')
      s,p,b = self.parseDetID(int(detID))
      if int(SiPM)>7: side='R'
      else: side='L'
      if gelsides[p]==side: return 1
      else: return 0

   def dist2BarEnd(self, MuFilter, detID, nSides, pred, dim):
      if nSides == 2 and dim=='x':
      	MuFilter.GetPosition(detID, self.A, self.B)
      	dxL=pred-B.x()
      	dxR=pred-A.x()
      	return dxL, dxR

      elif nSides==2 and dim=='y':
      	MuFilter.GetPosition(detID, self.A, self.B)
      	dyT=A.y()-pred
      	return dyT		

      elif nSides == 1 and dim=='y':
      	MuFilter.GetPosition(detID,self.A, self.B)
      	dyT=A.y()-pred
      	return dyT

      elif nSides == 1 and dim=='y':
      	MuFilter.GetPosition(detID, self.A, self.B)
      	dyT=A.y()-pred
      	return dyT   

   def parseDetID(self, detID):
      if not isinstance(detID, int): detID=int(detID)
      subsystem=detID//10000
      if subsystem in (1,2,3):
         plane=detID%10000//1000
         bar=detID%1000
         return subsystem, plane, bar

   def MakeDetID(self, fixed):
      fixed_subsystem, fixed_plane, fixed_bar = fixed
      if fixed_subsystem in (1,2,3):
         return int(f'{str(fixed_subsystem)}{str(fixed_plane)}{str(fixed_bar).zfill(3)}')

   def MakeFixedCh(self, fixed):
      fixed_subsystem, fixed_plane, fixed_bar, fixed_SiPM = fixed
      if fixed_subsystem in (1,2,3):
         return f'{str(fixed_subsystem)}{str(fixed_plane)}{str(fixed_bar).zfill(3)}_{str(fixed_SiPM)}'

   def GetDSHaverage(self, hits, nPlanes=3):
      stations={k:{} for k in range(7)}
      for i,hit in enumerate(hits):
         detID=hit.GetDetectorID()
         s,p,b=self.parseDetID(detID)
         if s!=3:continue
         if b>59: continue # The bar IDs in DS4V start at 60 so the vertical bars are not contaminating here.
         stations[p][i]=hit
      # if not all( len(stations[p*2])==1 for p in range(nPlanes) ): # Either the 2 or 3 horizontal DS planes.
      if not all( [len(stations[p])==1 for p in range(nPlanes)] ):
         return -999.
      # ts=[0.,0] # 1st element to sum tdcs, 2nd element to count channels (always equal to 6)
      total=0.
      for p in range(nPlanes):
         hit=stations[p][list(stations[p].keys())[0]]
         tdcs=hit.GetAllTimes()
         if len(tdcs) != 2: return -999. # Event selection: only accept events with BOTH SiPMs firing on 1 scintillator
         # for all 3 horizontal DS planes.
         for item in tdcs: 
            SiPM, tdc = item
            total+=tdc 
      dsh_average=total/6 # Due to event selection in L146, there are always 6 channels over which to be averaged.
      return dsh_average

   def ATLAStrack(self, hits, DST0):
      for hit in hits:
         detID=hit.GetDetectorID()
         s,p,b=self.parseDetID(detID)
         if s==2 and p==0:
            US1hit=hit
            break
      averageUS1TDC=self.GetAverageTDC(US1hit)
      if averageUS1TDC < DST0: return True 
      else: return False

   def fit_langau(self, hist,o,bmin,bmax):
      params = {0:'Width(scale)',1:'mostProbable',2:'norm',3:'sigma'}
      F = ROOT.TF1('langau',langaufun,0,200,4)
      for p in params: F.SetParName(p,params[p])
      rc = hist.Fit('landau','S'+o,'',bmin,bmax)
      res = rc.Get()
      if not res: return res
      F.SetParameter(2,res.Parameter(0))
      F.SetParameter(1,res.Parameter(1))
      F.SetParameter(0,res.Parameter(2))
      F.SetParameter(3,res.Parameter(2))
      F.SetParLimits(0,0,10)
      F.SetParLimits(1,0,100)
      F.SetParLimits(3,0,10)

      rc = hist.Fit(F,'S'+o,'',bmin,bmax)
      res = rc.Get()
      return res

   def langaufun(self, x,par):
      #Fit parameters:
      #par[0]=Width (scale) parameter of Landau density
      #par[1]=Most Probable (MP, location) parameter of Landau density
      #par[2]=Total area (integral -inf to inf, normalization constant)
      #par[3]=Width (sigma) of convoluted Gaussian function
      #
      #In the Landau distribution (represented by the CERNLIB approximation),
      #the maximum is located at x=-0.22278298 with the location parameter=0.
      #This shift is corrected within this function, so that the actual
      #maximum is identical to the MP parameter.
      #
      # Numeric constants
      invsq2pi = 0.3989422804014   # (2 pi)^(-1/2)
      mpshift  = -0.22278298       # Landau maximum location
      #
      # Control constants
      np = 100.0      # number of convolution steps
      sc =   5.0      # convolution extends to +-sc Gaussian sigmas
      #
      # Variables
      summe = 0.0
      #
      # MP shift correction
      mpc = par[1] - mpshift * par[0]
      #
      # Range of convolution integral
      xlow = x[0] - sc * par[3]
      xupp = x[0] + sc * par[3]
      #
      step = (xupp-xlow) / np
      #
      # Convolution integral of Landau and Gaussian by sum
      i=1.0
      if par[0]==0 or par[3]==0: return 9999
      while i<=np/2:
         i+=1
         xx = xlow + (i-.5) * step
         fland = ROOT.TMath.Landau(xx,mpc,par[0]) / par[0]
         summe += fland * ROOT.TMath.Gaus(x[0],xx,par[3])
         #
         xx = xupp - (i-.5) * step
         fland = ROOT.TMath.Landau(xx,mpc,par[0]) / par[0]
         summe += fland * ROOT.TMath.Gaus(x[0],xx,par[3])
         #
      return (par[2] * step * summe * invsq2pi / par[3])

   def GetAverageTDC(self, mufiHit, side='both'):
      value=[0, 0]
      count=[0, 0]
      nSiPMs=mufiHit.GetnSiPMs()
      times=mufiHit.GetAllTimes()
      s, p, b = self.parseDetID(mufiHit.GetDetectorID())
      for element in times:
         SiPM, time = element
         if s==2 and self.IsSmallSiPMchannel(SiPM):continue
         if SiPM<nSiPMs:
            value[0]+=time
            count[0]+=1
         else:
            value[1]+=time
            count[1]+=1
      if s == 2:
         if count[0] != 0 and count[1] != 0:
            if side=='both':
               average = 0.5*(value[0]/count[0]+value[1]/count[1])
               return average
            if side=='L':
               return value[0]/count[0]
            elif side == 'R':
               return value[1]/count[1]
         else: return -999.
      elif s == 3 and b < 60: 
         if side=='both':
            average = 0.5*(value[0]/count[0]+value[1]/count[1])
            return average
         elif side == 'L': return value[0]/count[0]
         elif side == 'R': return value[1]/count[1]

   def GetChannelVal(self, SiPM, chs):
      for entry in chs:
         fSiPM, val = entry
         if fSiPM == SiPM:
            return val
      return -999.

   def TotalChVal_sides(self, chs):
      left, right=0., 0.
      for entry in chs:
         SiPM, val = entry
         if SiPM<8: left+=val
         elif SiPM>=8: right+=val
      return left, right

   def GetOverallSiPMNumber(self, detID, SiPM):
      s, p, b = self.parseDetID(detID)
      if s==1: nSiPMs, SiPMs_plane=8, 56 # is it?
      elif s==2: nSiPMs, SiPMs_plane=16, 160
      elif s==3: nSiPMs, SiPMs_plane=1, 60 # wrong

      return SiPM+nSiPMs*b+p*SiPMs_plane   

   def OneHitPerSystem(self, hits, systems, Nfired=False):
   # def OneHitPerSystem(self, hits, systems, Nfired=False, verbose=False):
      verbose=self.verbose

      hitdict={}
      for s in systems: # systems always includes US and DS. Veto is included if a Scifi track is also formed.
         for p in range(self.systemAndPlanes[s]):
            key=10*s+p
            hitdict[key]=[]
      for i, hit in enumerate(hits):
         if not hit.isValid(): continue
         detID=hit.GetDetectorID()
         s, p, b = self.parseDetID(detID)
         if s not in systems: continue
         key=10*s+p
         hitdict[key].append(detID)

      #### For returning fraction of planes with 1 scintillator firing.
      if Nfired:
         counter=0
         totalplanes=12 if systems==(2,3) else 14
         
         for key in hitdict:
            if key<30 and len(hitdict[key])==1: counter+=1
            elif key>=30 and key<33: 
               if len(hitdict[key])!=2: continue
               bars=sorted([ self.parseDetID(hit)[2] for hit in hitdict[key] ]) # List of bar numbers sorted by number. 
               if bars[0]<60: counter+=1
               if bars[1]>59: counter+=1
            elif key==33 and len(hitdict[key])==1:counter+=1
         return counter/totalplanes 

      #### For requiring an event to have exactly 1 scintillator firing.
      for key in hitdict:
         hits=hitdict[key]
         if key<30 and len(hitdict[key]) != 1: return False # Veto and US

         elif key>=30 and key<33: # DS0 -> DS2
            if len(hits)!=2:
               if verbose: print(f'{key}: {hits}')
               return False 
            bars=sorted([self.parseDetID(hit)[2] for hit in hits]) # List of bar numbers sorted by number. 
            # [horizontal, vertical]
            if not bars[0]<60 and bars[1]>59: 
               if verbose: print(f'bars fired for plane {key}: {bars}')
               return False
         elif key==33: # DS3 only vertical bars
            if not len(hits)==1:
               return False
            # if not DSVcheck(hitdict[key]): return False
         return True

   def InAcceptance(self, pos, mom, subsystem, geoobject, zPos):

      limits=self.GetSubsystemZlimits(subsystem, zPos)
      if limits==0:return 0

      for z in limits:
         lam=(z-pos.z())/mom.z()
         Ex=ROOT.TVector3(pos.x()+lam*mom.x(), pos.y()+lam*mom.y(), pos.z()+lam*mom.z())
         xmin,xmax,ymin,ymax=self.GetSubsystemXYlimits(subsystem, geoobject)
         inacc=xmin<Ex.x()<xmax and ymin<Ex.y()<ymax
         # print(f'{xmin}<{Ex.x()}<{xmax}, {ymin}<{Ex.y()}<{ymax}, {inacc}')
         if not inacc: return False
      return True

   def GetSubsystemZlimits(self, subsystem, zPos):
      if subsystem==0:
         pass
      elif subsystem==1:
         zmin, zmax=zPos['MuFilter'][10], zPos['MuFilter'][11]
      elif subsystem==2:
         zmin, zmax=zPos['MuFilter'][20], zPos['MuFilter'][24]
      elif subsystem==3:
         zmin, zmax=zPos['MuFilter'][30], zPos['MuFilter'][36]
      else: return 0
      return zmin, zmax

   def GetSubsystemXYlimits(self, subsystem,geoobject):
      if subsystem==0:
         pass
      elif subsystem==1:
         geoobject.GetPosition(10000, A, B)
         xmin, xmax = B.x(), A.x()
         ymin = A.y() - geoobject.GetConfParF('MuFilter/VetoBarY')/2.
         geoobject.GetPosition(10006, A, B)
         ymax = A.y() + geoobject.GetConfParF('MuFilter/VetoBarY')/2.
      elif subsystem==2:
         geoobject.GetPosition(20000, A, B)
         xmin, xmax = B.x(), A.x()
         ymin = A.y() - geoobject.GetConfParF('MuFilter/UpstreamBarY')/2.
         geoobject.GetPosition(20009, A, B)
         ymax = A.y() + geoobject.GetConfParF('MuFilter/UpstreamBarY')/2.
      elif subsystem==3: 
         pass
      return xmin, xmax, ymin, ymax

   def getNUSPlanes(self, hits):

      res=0

      USPlanes={k:0 for k in range(5)}
      for i, hit in enumerate(hits):
         if not hit.isValid():continue
         detID=hit.GetDetectorID()
         s, p, b=self.parseDetID(detID)
         if s != 2: continue
         USPlanes[p]+=1
      for plane in range(5):
         if USPlanes[plane]==1: res+=1
      return res	

   def delta_min_t(self, aHit):
      times = aHit.GetAllTimes()
      if len(times)==0: return -999.
      nSiPMs = aHit.GetnSiPMs()
      ts_L, ts_R = [], []
      for channel in times:
         SiPM, time = channel
         if SiPM<nSiPMs: ts_L.append(time)
         elif SiPM>=nSiPMs: ts_R.append(time)  
      return min(ts_L)-min(ts_R)

   def GetnFiredSiPMs(self, aHit):
      nSiPMs=aHit.GetnSiPMs()

      nFiredSiPMs_left=0
      nFiredSiPMs_right=0
      channels=aHit.GetAllSignals()
      for ch in channels:
      	SiPM, qdc = ch 
      	if SiPM<nSiPMs: nFiredSiPMs_left+=1
      	elif SiPM>=nSiPMs: nFiredSiPMs_right+=1
      return nFiredSiPMs_left, nFiredSiPMs_right 

   def GetSiPMNumberInSystem_LandR(self, detID, SiPM): # 20000 SiPM 8 -> 8
      if not isinstance(SiPM, int): SiPM=int(SiPM)
      s, p, b = self.parseDetID(int(detID))
      if s==1:
         nSiPMs, SiPMs_plane=16, 112 # is it? 
         return int(SiPM)+nSiPMs*b+p*SiPMs_plane
      elif s==2:
         nSiPMs, SiPMs_plane=16, 160
         return SiPM+nSiPMs*b+p*SiPMs_plane

      elif s==3: # Count left and right horizontal SiPMs consecutively
         nSiPMs_bar_hor, nSiPMs_bar_ver=2, 1
         nSiPMs_plane_hor, nSiPMs_plane_ver=120, 60
         # p=verticalBarDict[p]
         # tmp=(p-1) if p!=0 else 0
         # if SiPM>59:
         #    nsipm=tmp*nSiPMs_plane_hor+p*nSiPMs_plane_ver+nSiPMs_ver*b+SiPM
         # elif SiPM<60:
         #    nsipm=p*nSiPMs_plane_ver+tmp*nSiPMs_plane_hor+nSiPMs_hor*b+SiPM
         # return nsipm
         if b>59:
            tmp=p+1
            horizontalSiPMs=(p+1)*nSiPMs_plane_hor
            verticalSiPMs=p*nSiPMs_plane_ver
            total=horizontalSiPMs+verticalSiPMs+(b-60)*nSiPMs_bar_ver
         elif b<60 and p==3:
            horizontalSiPMs=p*nSiPMs_plane_hor
            verticalSiPMs=p*nSiPMs_plane_ver
            total=horizontalSiPMs+verticalSiPMs+b*nSiPMs_bar_ver
         elif b<60 and p!=3:
            horizontalSiPMs=p*nSiPMs_plane_hor
            verticalSiPMs=p*nSiPMs_plane_ver
            total=horizontalSiPMs+verticalSiPMs+b*nSiPMs_bar_hor+SiPM
         return total
   
   def GetSiPMNumberInPlane_LandR(self, detID, SiPM):
      s, p, b = self.parseDetID(detID)
      if s == 1: return SiPM*b*12
      if s == 2:  return SiPM+b*16 

   def GetSiPMNumberInPlane_LTR(detID, SiPM):
      s, p, b = self.parseDetID(detID)
      if s != 2: 
         print('AAAAAAAHHHHHH')
      return SiPM+b*8

   def GetSiPMNumberInSystem_LTR(self, detID, SiPM): # 20000 SiPM 8 -> 400
      if not isinstance(SiPM, int): SiPM=int(SiPM)
      s, p, b = self.parseDetID(detID)
      if s==1: nSiPMs, SiPMs_plane=6, 52 # is it?
      elif s==2: nSiPMs, SiPMs_plane=8, 80
      elif s==3: nSiPMs, SiPMs_plane=1, 60 # wrong 

      if SiPM<nSiPMs:
         return SiPM+nSiPMs*b+p*SiPMs_plane
      elif SiPM>=nSiPMs:
         SiPM_r=400+SiPM%nSiPMs
         return SiPM_r+nSiPMs*b+p*SiPMs_plane

   def GetDeltaT(self, times, one_channel=None):
      # nSiPMs=aHit.GetnSiPMs()
      nSiPMs=8
      mean = [0,0]
      count = [0,0]
      # channels = aHit.GetAllTimes()
      for ch in times:
         SiPM, val = ch
         if one_channel != None:
            if not (SiPM == one_channel or SiPM == one_channel+nSiPMs): continue
         if self.IsSmallSiPMchannel(SiPM): continue
         if SiPM < nSiPMs:
            mean[0]+=val
            count[0]+=1
         else:
            mean[1]+=val
            count[1]+=1
      # print(count, mean)
      if count[0] != 0 and count[1] != 0:
         return (mean[0]/count[0]-mean[1]/count[1])/2.
      else: return -999.

   def GetdtCalc(self, xpred, L, cs):
      left, right=list(filter(lambda x : x[0]<8, cs)), list(filter(lambda x : x[0]>7, cs))
      NL, NR=len(left), len(right)
      if NL == 0 or NR == 0: return -999.
      sumOfInverses=lambda x : sum( [1/i[1] for i in x] )
      return xpred/NL*sumOfInverses(left) - (L-xpred)/NR*sumOfInverses(right)

   def Getcscint(self, runNr, fixed_ch, state):

      iteration=0 if state=='uncorrected' else 1
      if not os.path.exists(f'{self.path}cscintvalues/run{runNr}/cscint_{fixed_ch}.csv'): 
         return -999
      with open(f'{self.path}cscintvalues/run{runNr}/cscint_{fixed_ch}.csv', 'r') as handle:
         reader=csv.reader(handle)
         alldata=[row for row in reader]
         if len(alldata)<iteration+1 or len(alldata) == 0: 
            print(f'{fixed_ch} Len issue')
            return -999
         try:
            data=alldata[iteration]
         except IndexError:
            print(f'{fixed_ch} IndexError')
         return (float(data[1]), float(data[2]))

   def Getcscint_chi2pNDF_info(self, runNr,fixed_ch,state):
      iteration=0 if state=='uncorrected' else 1
      with open(f'{self.path}cscintvalues/run{runNr}/cscint_{fixed_ch}.csv', 'r') as handle:
         reader=csv.reader(handle)
         alldata=[row for row in reader]
         if len(alldata)<iteration+1: return -999.
         try:
            data=alldata[iteration]
         except IndexError:
            print(f'{fixed_ch} IndexError')
      return float(data[-2]),int(data[-1])

   def Getcscint_chi2pNDF(self, runNr,fixed_ch,state):
      iteration=0 if state=='uncorrected' else 1
      if not os.path.exists(f'{self.path}cscintvalues/run{runNr}/cscint_{fixed_ch}.csv'): return -999.
      with open(f'{self.path}cscintvalues/run{runNr}/cscint_{fixed_ch}.csv', 'r') as handle:
         reader=csv.reader(handle)
         alldata=[row for row in reader]
         if len(alldata)<iteration+1: return -999.
         try:
            data=alldata[iteration]
         except IndexError:
            print(f'{fixed_ch} IndexError')
      return float(data[-2])/int(data[-1])

   def Makecscintdict(self, runNr, subsystem, state):
      iteration=0 if state=='uncorrected' else 1
      path=f'{self.path}/cscintvalues/run{runNr}/'
      res={}
      for filename in os.listdir(path):
         #fixed_ch=filename[filename.find(str(subsystem)):filename.find('.csv')]
         fixed_ch=f'{filename.split("_")[1]}_{filename.split("_")[2].split(".")[0]}'
         detID=int(fixed_ch.split('_')[0])
         s,p,b = self.parseDetID(detID)
         if s!=subsystem: continue
         if s==3 and DSVcheck(detID): 
            print(detID)
            continue
         with open(path+filename, 'r') as handle:
            reader=csv.reader(handle)
            alldata=[row for row in reader]
            if len(alldata)<iteration or alldata==[]:
               print(filename)
               continue
            data=alldata[iteration]
         res[fixed_ch]=float(data[1])
      sorted_tuples=sorted(res.items(), key=lambda x:x[1])
      sorted_d={k:v for k,v in sorted_tuples}
      return sorted_d

   def Maketimeresolutiondict(self, runNr, subsystem, state):
      iteration=0 if state=='uncorrected' else 1
      path=f'{self.path}TimeResolution/run{runNr}/'
      res={}
      for filename in os.listdir(path):
         fixed_ch=f'{filename.split("_")[1]}_{filename.split("_")[2].split(".")[0]}'
         if fixed_ch[0]!=str(subsystem): continue
         with open(path+filename, 'r') as handle:
            reader=csv.reader(handle)
            alldata=[row for row in reader]
            if len(alldata)<iteration or alldata==[]:
               print(filename)
               continue
            data=alldata[iteration]
         res[fixed_ch]=(float(data[1]), float(data[2]))
      sorted_tuples=sorted(res.items(), key=lambda x:x[1][0])
      sorted_d={k:v for k,v in sorted_tuples}
      return sorted_d   

   def GetPolyParams(self, runNr, fixed_ch, n, state='uncorrected'):
      iteration=0 if state=='uncorrected' else 1
      filelengths={1:11, 2:13, 3:15, 4:13, 5:9}
      if not os.path.exists(f'{self.path}Polyparams/run{runNr}/polyparams{n}_{fixed_ch}.csv'): return -999.
      with open(f'{self.path}Polyparams/run{runNr}/polyparams{n}_{fixed_ch}.csv', 'r') as f:
         reader=csv.reader(f)
         alldata=[r for r in reader]
         if len(alldata)==0:return -999
         data=alldata[iteration]
      if len(data)!=filelengths[n]: return -999
        
      params=[float(i) for i in data[1:-2]]
      limits=[float(i) for i in data[-2:]]
      return params,limits

   def Gettimeresolution(self, runNr, fixed_ch, state):
      iteration=0 if state=='uncorrected' else 1
      fname=f'{self.path}TimeResolution/run{runNr}/timeresolution_{fixed_ch}.csv'
      if not os.path.exists(fname): return -999
      with open(fname, 'r') as f:
         reader=csv.reader(f)
         alldata=[row for row in reader]
         if len(alldata)<iteration+1: return -999.
         data=alldata[iteration]
         # if math.isnan(float(data[0])): return -999.
      timeresolution=float(data[1]), float(data[2])
      return timeresolution

   def FitForMPV(self, runNr, fixed_ch, state):
      iteration=0 if state=='uncorrected' else 1
      fname=f'{self.path}rootfiles/run{runNr}/timewalk_{fixed_ch}.root'
      if not os.path.exists(fname): return -999.
      f=ROOT.TFile.Open(fname, 'READ')
      histname=f'dtvqdc_{fixed_ch}_{state}'
      if not histname in [k.GetName() for k in f.GetListOfKeys()]: return -999
      tmp=f.Get(histname)
      hist=tmp.Clone()
      hist.SetDirectory(ROOT.gROOT)
      f.Close()
      xproj=hist.ProjectionX()
      xmin,xmax=xproj.GetXaxis().GetXmin(), xproj.GetXaxis().GetXmax()
      res=fit_langau(xproj, 'S Q', xmin, xmax)
      MPV=res.Parameter(1)
      MPV_err=res.ParError(1)
      chi2, NDF= res.Chi2(), res.Ndf()
      return (MPV, MPV_err, chi2, NDF)

   def Getchi2_info(self, runNr, fixed_ch, n, state):
      iteration=0 if state=='uncorrected' else 1

      fname=f'{self.path}chi2s/run{runNr}/chi2s{n}_{fixed_ch}.csv'
      if not os.path.exists(fname): return -999.
      with open(fname, 'r') as handle:
         reader=csv.reader(handle)
         alldata=[row for row in reader]
         if len(alldata)==0 or len(alldata)<iteration: return -999.
         data=alldata[iteration-1]
         try:
            chi2info=str(data[0]), float(data[1]), int(data[2])
         except ValueError:
            print(f'Non-integer NDF for {fixed_ch}')
            return -999.
      return chi2info

   def Getchi2pNDF(self, runNr, fixed_ch, n, state):

      iteration=0 if state=='uncorrected' else 1   
      fname=f'{self.path}chi2s/run{runNr}/chi2s{n}_{fixed_ch}.csv'
      if not os.path.exists(fname): return -999.
      with open(fname, 'r') as handle:
         reader=csv.reader(handle)
         alldata=[row for row in reader]
         if len(alldata)==0 or len(alldata)<iteration: return -999.
         data=alldata[iteration-1]
         try:
            chi2pNDF=float(data[1])/int(data[2])
         except ValueError:
            print(f'Non-integer NDF for {fixed_ch}')
            return -999.
      return chi2pNDF

   def GetNDF(self, runNr, fixed_ch, iteration):
      fname=self.path+'chi2s/run'+str(runNr)+'/chi2s_'+fixed_ch+'.csv'
      if not os.path.exists(fname): return -999.
      with open(fname, 'r') as handle:
         reader=csv.reader(handle)
         alldata=[row for row in reader]
         if len(alldata)==0 or len(alldata)<iteration: return -999.
         data=alldata[iteration-1]
      return int(data[2])

   def GetBadchi2pNDFdict(self, runNr, subsystem, state):

      iteration=0 if state=='uncorrected' else 1
      path=f'{self.path}/chi2s/run{runNr}/'
      res={}
      for filename in os.listdir(path):
         #fixed_ch=filename[filename.find(str(subsystem)):filename.find('.csv')]
         fixed_ch=f'{filename.split("_")[1]}_{filename.split("_")[2].split(".")[0]}'
         if fixed_ch[0]!=str(subsystem): continue
         with open(path+filename, 'r') as handle:
            reader=csv.reader(handle)
            alldata=[row for row in reader]
            if len(alldata)<iteration:
               print(filename)
               continue
            data=alldata[iteration-1]
            #chi2pNDF=data.pop()
            chi2pNDF=float(data[1])/int(data[2])
         res[fixed_ch]=chi2pNDF
      sorted_tuples=sorted(res.items(), key=lambda x:x[1])
      sorted_d={k:v for k,v in sorted_tuples}
      return sorted_d

   def Makechi2pNDFdict(self, runNr, subsystem, n, state):

      iteration=0 if state=='uncorrected' else 1
      path=f'{self.path}/chi2s/run{runNr}/'
      res={}
      for filename in os.listdir(path):
         if filename.split('_')[0].find(f'chi2s{n}') ==-1: continue
         #fixed_ch=filename[filename.find(str(subsystem)):filename.find('.csv')]
         fixed_ch=f'{filename.split("_")[1]}_{filename.split("_")[2].split(".")[0]}'
         if fixed_ch[0]!=str(subsystem): continue
         with open(path+filename, 'r') as handle:
            reader=csv.reader(handle)
            alldata=[row for row in reader]
            if len(alldata)==0 or len(alldata)<iteration: return -999.
            if len(alldata)<iteration:
               print(filename)
               continue
            data=alldata[iteration-1]
            chi2pNDF=float(data[1])/int(data[2])
         res[fixed_ch]=chi2pNDF
      sorted_tuples=sorted(res.items(), key=lambda x:x[1])
      sorted_d={k:v for k,v in sorted_tuples}
      return sorted_d

   def GetAttenuationLengthData(self, runNr, fixed_ch):
      fname=self.path+'attenuationlengths/run'+str(runNr)+'/csvfiles/attenuationlength_'+fixed_ch+'.csv'
      if not os.path.exists(fname): return -999.
      with open(fname, 'r') as handle:
         reader=csv.reader(handle)
         alldata=[row for row in reader]
         data=alldata[0]
         return data

   def GetXcalculated(self, dt, L, cs, wanted=None):

      left, right=list(filter(lambda x : x[0]<8, cs)), list(filter(lambda x : x[0]>7, cs))
      NL, NR=len(left), len(right)
      if NL==0 or NR == 0: return -999.
      sumOfInverses=lambda x : sum( [1/i[1] for i in x] )
      A, B = 1/NL*sumOfInverses(left), 1/NR*sumOfInverses(right)
      xcalc = (dt+L*B)/(A+B)

      return xcalc

   def GetMPV(self, runNr, fixed_ch, iteration):
      fname=f'{self.path}MPVs/run{runNr}/MPV_{fixed_ch}.csv'
      if not os.path.exists(fname): return -999.
      with open(fname, 'r') as h:
         reader=csv.reader(h)
         data=[row for row in reader]
         if data==[]: return -999. # To be investigated
         res=data[iteration-1]
      return float(res[0])

   def GetToFcorrection(self, SiPM, xpred, cs, xref):
      c_SiPM=float(cs[0])
      ToFcorrection=abs((xpred-xref)/c_SiPM)
      return ToFcorrection

   def ApplyToFCorrection(self, SiPM, ToF, clock, xpred, xref):

      time=clock*6.25
      if SiPM<8:
         if xpred >= xref: corrected_t = time-ToF
         else: corrected_t = time+ToF
      else: 
         if xpred >= xref: corrected_t = time+ToF
         else: corrected_t = time-ToF
      return (SiPM, corrected_t)
   
   def correct_ToF(self, SiPM, clock, xpred, cs, xref):

      # fixed_subsystem, fixed_plane, fixed_bar, fixed_SiPM = fixed
      time=clock*6.25
      c_SiPM=float(cs[0])
      ToFcorrection=abs((xpred-xref)/c_SiPM)
      # Does not work for DS! 
      if SiPM<8:
         if xpred >= xref: corrected_t = time - ToFcorrection
         else: corrected_t = time + ToFcorrection
      else: 
         if xpred >= xref: corrected_t = time + ToFcorrection
         else: corrected_t = time - ToFcorrection
      return (SiPM, corrected_t)

   def GetCutDistributions(self, runNr, distmodes=('dy', 'slopes', 'nSiPMs'), nStations=2):
      Allmodes=('dy', 'nSiPMs', 'slopes')
      filename=f'{self.path}rootfiles/run{runNr}/SelectionCriteria.root'
      if not os.path.exists(filename): filename=f'{self.path}rootfiles/run005097/SelectionCriteria.root'

      if isinstance(distmodes, str):
         distmodes=(distmodes,)

      for distmode in distmodes:
         if distmode not in Allmodes:
            print('Use a valid mode.')
            return 0

      f=ROOT.TFile.Open(filename, 'READ')
      for distmode in distmodes:
         if not hasattr(f, distmode):
            f.Close()
            filename=f'{self.path}rootfiles/run005097/SelectionCriteria.root'
            f=ROOT.TFile.Open(filename, 'READ')
            break
      dists={}

      for distmode in distmodes:
         if distmode=='dy' or distmode=='nSiPMs':
            for s in (1,2):
               for p in range(self.systemAndPlanes[s]):
                  key=str(s*10+p)
                  name=f'{distmode}_{key}_{nStations}stations'
                  if not hasattr(f, name):
                     print(f'No hist {name}')
                     f.Close()
                     return -999
                  hist=f.Get(name).Clone()
                  hist.SetDirectory(ROOT.gROOT)
                  dists[name]=hist

         elif distmode=='slopes':
            if not hasattr(f, name):
               f.Close()
               print(f'No hist {name}')
               return -999
            hist=f.Get(f'{distmode}_{nStations}stations').Clone()
            hist.SetDirectory(ROOT.gROOT)
            dists[distmode]=hist

      f.Close()
      return dists

   def GetEntriesInHist(self, runNr, fixed_ch, mode, state):
       f=ROOT.TFile.Open(f'{self.path}/rootfiles/run{runNr}/timewalk_{fixed_ch}.root', 'READ')
       hist=f.Get(f'{mode}_{fixed_ch}_{state}')
       # hist.SetDirectory(ROOT.gROOT)
       entries=hist.GetEntries()
       f.Close()
       return entries

   def GetPolyParamRanges(self, runNr, n, subsystem, iteration):
      params={i:[0.,0.] for i in ('A', 'B', 'C')}
      path=f'{self.path}Polyparams/run{runNr}/'
      files=[i for i in os.listdir(path) if int(i.split('_')[1][0])==subsystem]
      counter=0
      for f in files:
         if f.split('_')[0].find(str(n)) == -1:continue
         #if idx>1: break
         fixed_ch=f"{f.split('_')[1]}_{f.split('.')[0].split('_')[-1]}"
         tmp = self.GetPolyParams(runNr, fixed_ch, n, iteration)
         if isinstance(tmp, int): continue
         else: fps=tmp[0]
         if fps==-999.: 
            print(f'{fixed_ch} has no poly params')
            continue

         if n==1:
            vals={p:fps[i*2] for i,p in enumerate(('A', 'B', 'C'))}
            if counter==0:
                for i in ('A', 'B', 'C'): params[i]=[vals[i], 1.1*vals[i]]
                # for idx,x in enumerate(('A', 'B', 'C')): params[x]=[fps[idx*2], fps[idx*2+1]]
            else: 
               for i in ('A', 'B', 'C'):
                  if vals[i] < params[i][0]: params[i][0]=vals[i]
                  if vals[i] > params[i][1]: params[i][1]=1.1*vals[i]
            if n==5:
               vals={p:fps[i*2] for i,p in enumerate(('A', 'B', 'C'))}
               if counter==0:
                   for i in ('A', 'B', 'C'): params[i]=[vals[i], 1.1*vals[i]]
                   # for idx,x in enumerate(('A', 'B', 'C')): params[x]=[fps[idx*2], fps[idx*2+1]]
               else: 
                  for i in ('A', 'B', 'C'):
                     if vals[i] < params[i][0]: params[i][0]=vals[i]
                     if vals[i] > params[i][1]: params[i][1]=1.1*vals[i]         
            counter+=1

      return params

   def MakePolyParamDicts(self, runNr, subsystem, iteration):
      #A,B,C={}, {}, {}
      params={i:[] for i in ('A', 'B', 'C')}
      path=f'{self.path}Polyparams/run{runNr}/'
      files=[i for i in os.listdir(path) if int(i.split('_')[1][0])==subsystem]
      for idx, f in enumerate(files):
         fixed_ch=f"{f.split('_')[1]}_{f.split('.')[0].split('_')[-1]}"
         fps=self.GetPolyParams(runNr, fixed_ch, iteration)
         if fps==-999.: continue
         vals=[fps[i*2] for i in range(3)]
         [params[p].append(vals[i]) for i,p in enumerate(params)]
         ranges={i:[0,0] for i in ('A', 'B', 'C')}
      return params

   def GetCanvas(self, runNr, fixed_ch, mode, iteration):
      filename=f'{self.path}rootfiles/run{runNr}/timewalk_{fixed_ch}.root'
      infile=ROOT.TFile.Open(filename, 'read')
      name=f'{mode}_{fixed_ch}_{iteration}'
      if not hasattr(infile, name):
         print(f'No canvas with name: {name} available in file')
         return 
      og_canv=infile.Get(name)
      canv.og_canv.Clone()
      canv.SetDirectory(ROOT.gROOT)
      infile.Close()
      return canv

   def WriteCanvasesToFile(self, runNr, mode, iteration, subsystems=(1,2)):
      pathtodir=f'{self.path}Results/{runNr}'
      if not os.path.exists(pathtodir):
         os.mkdir(pathtodir)
      outfilename=pathtodir+f'run{runNr}_results.root'
      outfile=ROOT.TFile.Open(outfilename, 'recreate')
      files = [f'{self.path}run{runNr}/timewalk_{self.MakeFixedCh((s,p,b,SiPM))}.root' for s in subsystems for p in range(self.systemAndPlanes[s]) for p in range(self.systemAndBars[s]) for b in range(self.systemAndBars[s]) for SiPM in self.systemAndSiPMs[s] ]
      for i, infile in enumerate(files):
         fixed_ch=infile[infile.find('timewalk')+len('timewalk_'):infile.find('.root')]
         canv=GetCanvas(runNr, fixed_ch, mode, iteration)
         if not canv: 
            print(f'No {mode} canvas for {fixed_ch}')
            continue
         outfile.WriteObject(canv, canv.GetName(), 'kOverwrite')
      print(f'{i+1} canvases written to outfilename')
      outfile.Close()


