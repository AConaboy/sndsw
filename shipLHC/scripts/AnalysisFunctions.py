import ROOT,os,csv, math

# class mufi_analysis: # I could make several classes for generic functions like parseDetID and more specific ones for 
	# mapping SiPM channels to the whole system
afswork='/afs/cern.ch/work/a/aconsnd/Timing/'
afsuser='/afs/cern.ch/user/a/aconsnd/twfiles/'
A, B=ROOT.TVector3(), ROOT.TVector3()
systemAndPlanes = {1:2,2:5,3:6}
verticalBarDict={0:1, 1:3, 2:5, 3:6}
gelsides={0:'R', 1:'L', 2:'R', 3:'L', 4:'L'}

def BuildBarLengths(MuFilter):
   Vetobarlength = MuFilter.GetConfParF('MuFilter/VetoBarX')
   USbarlength = MuFilter.GetConfParF('MuFilter/UpstreamBarX')
   DSbarlength_hor = MuFilter.GetConfParF('MuFilter/DownstreamBarX')
   DSbarlength_vert = MuFilter.GetConfParF('MuFilter/UpstreamBarY_ver')
   barlengths={1:Vetobarlength, 2:USbarlength, 3:DSbarlength_hor, 4:DSbarlength_vert}
   return barlengths

def BuildzPos(MuFilter):
   zPos={}
   for s in systemAndPlanes:
      for plane in range(systemAndPlanes[s]):
         bar = 4
         p = plane
         if s==3 and plane%2==0:  
            bar = 90
            p = plane//2
         if s==3 and plane%2==1:
            bar = 30
            p = plane//2
         MuFilter.GetPosition(s*10000+p*1000+bar,A,B)
         zPos[s*10+plane] = (A.Z()+B.Z())/2.	
   return zPos

def IsSmallSiPMchannel(i):
	if i==2 or i==5 or i==10 or i==13: return True

	else: return False

def IsGel(fixed_ch):
   detID, SiPM = fixed_ch.split('_')
   s,p,b = parseDetID(int(detID))
   if int(SiPM)>7: side='R'
   else: side='L'
   if gelsides[p]==side: return 1
   else: return 0

def dist2BarEnd(MuFilter, detID, nSides, pred, dim):
	if nSides == 2 and dim=='x':
		MuFilter.GetPosition(detID, A, B)
		dxL=pred-B.x()
		dxR=pred-A.x()
		return dxL, dxR

	elif nSides==2 and dim=='y':
		MuFilter.GetPosition(detID, A, B)
		dyT=A.y()-pred
		return dyT		

	elif nSides == 1 and dim=='y':
		MuFilter.GetPosition(detID, A, B)
		dyT=A.y()-pred
		return dyT

	elif nSides == 1 and dim=='y':
		MuFilter.GetPosition(detID, A, B)
		dyT=A.y()-pred
		return dyT   

def parseDetID(detID):
   subsystem=detID//10000
   if subsystem ==1 or subsystem==2:
      plane=detID%10000//1000
      bar=detID%1000
      return subsystem, plane, bar
   if subsystem == 3:
      bar=detID%1000
      if bar>59:
         plane=verticalBarDict[detID%10000//1000]
      elif bar<60:
         plane=2*detID%10000//1000
      return subsystem, plane, bar

def fit_langau(hist,o,bmin,bmax):
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

   rc = hist.Fit(F,'S','',bmin,bmax)
   res = rc.Get()
   return res

def langaufun(x,par):
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

	      
def GetAvgT(mufiHit, side='both'):
   value=[0, 0]
   count=[0, 0]
   nSiPMs=mufiHit.GetnSiPMs()
   times=mufiHit.GetAllTimes()
   s, p, b = parseDetID(mufiHit.GetDetectorID())
   for element in times:
      SiPM, time = element 
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

   # elif s == 3 and b > 59: 
   #    return

def GetChannelVal(SiPM, chs):
	for entry in chs:
	   fSiPM, val = entry
	   if fSiPM == SiPM:
	      return val
	return -999.

def TotalChVal_sides(chs):
   left, right=0., 0.
   for entry in chs:
      SiPM, val = entry
      if SiPM<8: left+=val
      elif SiPM>=8: right+=val
   return left, right

def GetOverallSiPMNumber(detID, SiPM):
   s, p, b = parseDetID(detID)
   if s==1: nSiPMs, SiPMs_plane=6, 52 # is it?
   elif s==2: nSiPMs, SiPMs_plane=16, 160
   elif s==3: nSiPMs, SiPMs_plane=1, 60 # wrong

   return SiPM+nSiPMs*b+p*SiPMs_plane   

def OneHitPerSystem(hits,subsystems):
    Planes={k:0 for k in range(5)}
    if 1 in subsystems: Planes.update({10+k:0 for k in range(2)})
    for i, hit in enumerate(hits):
        if not hit.isValid(): continue
        detID=hit.GetDetectorID() 
        subsystem, plane, bar = parseDetID(detID)
        if subsystem != 2: continue
        Planes[plane]+=1
      
    for plane in range(5):
        if Planes[plane] != 1:    return False
    return True

def getNUSPlanes(hits):
   
   res=0

   USPlanes={k:0 for k in range(5)}
   for i, hit in enumerate(hits):
      if not hit.isValid():continue
      detID=hit.GetDetectorID()
      s, p, b=parseDetID(detID)
      if s != 2: continue
      USPlanes[p]+=1
   for plane in range(5):
      if USPlanes[plane]==1: res+=1
   return res

def DS_track(DigiHits):
	# check for low occupancy and enough hits in DS stations
   stations = {}
   # for s in systemAndPlanes:
   for plane in range(systemAndPlanes[3]): 

      stations[30+plane] = {}
    # k=-1
   # for i, aHit in enumerate(eventTree.Digi_MuFilterHit):
   for i, aHit in enumerate(DigiHits):
      # k+=1
      if not aHit.isValid(): continue
      detID=aHit.GetDetectorID()
      subsystem, plane, bar = parseDetID(detID)
      if subsystem != 3: continue
      key=subsystem*10+plane
      stations[key][i]=aHit
   if not len(stations[30])*len(stations[31])*len(stations[32])*len(stations[33]) == 1: return -1 # If not 1 hit in each DS plane
	#	build trackCandidate
   hitlist = {}
   for p in range(30,34):
      k = list(stations[p].keys())[0]
      hitlist[k] = stations[p][k]
   theTrack = trackTask.fitTrack(hitlist)
   return theTrack	

def delta_min_t(aHit):
   times = aHit.GetAllTimes()
   if len(times)==0: return -999.
   nSiPMs = aHit.GetnSiPMs()
   ts_L, ts_R = [], []
   for channel in times:
      SiPM, time = channel
      if SiPM<nSiPMs: ts_L.append(time)
      elif SiPM>=nSiPMs: ts_R.append(time)  
   return min(ts_L)-min(ts_R)

def GetnFiredSiPMs(aHit):
	nSiPMs=aHit.GetnSiPMs()

	nFiredSiPMs_left=0
	nFiredSiPMs_right=0
	channels=aHit.GetAllSignals()
	for ch in channels:
		SiPM, qdc = ch 
		if SiPM<nSiPMs: nFiredSiPMs_left+=1
		elif SiPM>=nSiPMs: nFiredSiPMs_right+=1
	return nFiredSiPMs_left, nFiredSiPMs_right 

def GetBarSlice(L, sliceL, xpred):
   nslices=int(L/sliceL)
   slice_num = int(xpred//sliceL)
   return slice_num   

def GetSiPMNumberInSystem_LandR(detID, SiPM): # 20000 SiPM 8 -> 8
   s, p, b = parseDetID(detID)
   if s==1: nSiPMs, SiPMs_plane=6, 52 # is it?
   elif s==2:
      nSiPMs, SiPMs_plane=16, 160
      return SiPM+nSiPMs*b+p*SiPMs_plane
   elif s==3: # Count left and right horizontal SiPMs consecutively
      nSiPMs, SiPMs_hor_plane, SiPMs_vert_plane=1, 120, 60
      if p not in verticalPlanes:
         total_SiPM = (p//2)*SiPMs_hor_plane+(p//2)*SiPMs_vert_plane+SiPM+2*b
         return total_SiPM
      else:
         if p==1: return SiPMs_hor_plane+b 
         elif p==3: return 2*SiPMs_hor_plane+SiPMs_vert_plane+b
         elif p==5: return 3*SiPMs_hor_plane+2*SiPMs_vert_plane+b
         elif p==6: return 3*SiPMs_hor_plane+3*SiPMs_vert_plane+b    
         
def GetSiPMNumberInPlane_LandR(detID, SiPM):
   s, p, b = parseDetID(detID)
   if s != 2: 
      print('AAAAAAAHHHHHH')
   return SiPM+b*16 

def GetSiPMNumberInPlane_LTR(detID, SiPM):
   s, p, b = parseDetID(detID)
   if s != 2: 
      print('AAAAAAAHHHHHH')
   return SiPM+b*8

def GetSiPMNumberInSystem_LTR(detID, SiPM): # 20000 SiPM 8 -> 400
   s, p, b = parseDetID(detID)
   if s==1: nSiPMs, SiPMs_plane=6, 52 # is it?
   elif s==2: nSiPMs, SiPMs_plane=8, 80
   elif s==3: nSiPMs, SiPMs_plane=1, 60 # wrong 

   if SiPM<nSiPMs:
      return SiPM+nSiPMs*b+p*SiPMs_plane
   elif SiPM>=nSiPMs:
      SiPM_r=400+SiPM%nSiPMs
      return SiPM_r+nSiPMs*b+p*SiPMs_plane

def GetDSH_average(hits, nPlanes=2):
   stations={k:{} for k in range(7)}
   for i,hit in enumerate(hits):
      detID=hit.GetDetectorID()
      s,p,b=parseDetID(detID)
      if b>59: continue
      if s!=3:continue
      stations[p][i]=hit
   if not all( len(stations[p*2])==1 for p in range(nPlanes) ): # Either the 2 or 3 horizontal DS planes.
      return -999.
   ts=[0.,0]
   for i in range(nPlanes):
      p=i*2
      hit=stations[p][list(stations[p].keys())[0]]
      tdcs=hit.GetAllTimes()
      if len(tdcs) != 2: 
         return -999.
      for item in tdcs: 
         SiPM, tdc = item
         ts[0]+=tdc
         ts[1]+=1
   dsh_average=ts[0]/ts[1]
   return dsh_average

def GetDeltaT(times, one_channel=None):
   # nSiPMs=aHit.GetnSiPMs()
   nSiPMs=8
   mean = [0,0]
   count = [0,0]
   # channels = aHit.GetAllTimes()
   for ch in times:
      SiPM, val = ch
      if one_channel != None:
         if not (SiPM == one_channel or SiPM == one_channel+nSiPMs): continue
      if IsSmallSiPMchannel(SiPM): continue
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

def GetdtCalc(xpred, L, cs):
   left, right=list(filter(lambda x : x[0]<8, cs)), list(filter(lambda x : x[0]>7, cs))
   NL, NR=len(left), len(right)
   if NL == 0 or NR == 0: return -999.
   sumOfInverses=lambda x : sum( [1/i[1] for i in x] )
   return xpred/NL*sumOfInverses(left) - (L-xpred)/NR*sumOfInverses(right)

def Getcscint_i(runNr, fixed_ch, iteration):
   with open(afswork+'cscintvalues/run'+str(runNr)+'/cscint_'+fixed_ch+'.csv', 'r') as handle:
      reader=csv.reader(handle)
      alldata=[row for row in reader]
      if len(alldata)<iteration or len(alldata) == 0: return -999.
      else:
         data=alldata[iteration]
         return (float(data[1]), float(data[2]))

def Getcscint_chi2pNDF(runNr,fixed_ch, iteration):
   # with open(afswork+'rootfiles/run'+str(runNr)+'/cscintvalues/cscint_'+fixed_ch+'.csv', 'r') as handle:
   with open(afswork+'cscintvalues/run'+str(runNr)+'/cscint_'+fixed_ch+'.csv', 'r') as handle:
      reader=csv.reader(handle)
      alldata=[row for row in reader]
      if len(alldata)<iteration or len(alldata) == 0: return -999.
      else:
         data=alldata[iteration]
         return data[-1]

def GetAvgcscint_i(path, detID, SiPMs, iteration):
   cs=[]
   for SiPM in SiPMs:
      fixed_ch=('_').join( (str(detID), str(SiPM)) ) 
      val=Getcscint_i(path, fixed_ch, iteration)
      cs.append(val)
   
   cbar = sum([i[0] for i in cs])/len(cs)
   cbar_err = ROOT.TMath.Sqrt( sum( [i[1]**2 for i in cs] ) )
   return cbar, cbar_err

def GetLogParams_i(runNr, fixed_ch, iteration):
   if iteration==0:
      print('No log params calculated for iteration 0.')
      return 0
   with open(afswork+'Logparams/run'+str(runNr)+'/logparams_'+fixed_ch+'.csv', 'r') as handle:
      reader=csv.reader(handle)
      alldata=[row for row in reader]
      if len(alldata)<iteration-1 or len(alldata) == 0: return -999.
      else:
         data=alldata[iteration-1]
         return (float(data[1]),float(data[2]),float(data[3]),float(data[4]),float(data[5]),float(data[6]))

def GetPolyparams_i(runNr, fixed_ch, iteration):
   # print(iteration)
   if iteration==0:
      print('No log params calculated for iteration 0.')
      return 0
   fname=afswork+'Polyparams/run'+str(runNr)+'/polyparams_'+fixed_ch+'.csv'
   if not os.path.exists(fname): return -999.
   with open(fname, 'r') as handle:
      reader=csv.reader(handle)
      alldata=[row for row in reader]
      if len(alldata)<iteration-1 or len(alldata)==0: return -999.
      else:
         if iteration>len(alldata):
            iteration-=1
         data=alldata[iteration-1][1:]
         data.pop() # We don't need to get the chi2 here.
         return tuple( [float(i) for i in data] )

def Gettimeresolution(runNr, fixed_ch, iteration):
   appendixdict={0:'_uncorrected', 1:'_corrected'}
   fname=afswork+'TimeResolution/run'+str(runNr)+'/timeresolution_'+fixed_ch+appendixdict[iteration]+'.csv'
   if not os.path.exists(fname): return -999
   with open(fname, 'r') as f:
      reader=csv.reader(f)
      alldata=[row for row in reader]
      data=alldata[0]
      if math.isnan(float(data[0])): return -999.
   timeresolution=float(data[1]), float(data[2])
   return timeresolution

def FitForMPV(runNr, fixed_ch, iteration):
   
   fname=afswork+'rootfiles/run'+str(runNr)+'/timewalk_dists_'+fixed_ch+'_poly.root'
   if not os.path.exists(fname): return -999.
   f=ROOT.TFile.Open(fname, 'READ')

   histname='dtcorrectedvqdc_'+fixed_ch+'_iteration'+str(iteration)
   if not histname in [k.GetName() for k in f.GetListOfKeys()]: return -999
   hist=f.Get(histname)
   xproj=hist.ProjectionX()
   #maxbin=xproj. # Currently all QDC histogram ranges are 0, 60 but this needs to change!
   res=fit_langau(xproj, 'Q', 0, 60)
   MPV=res.Parameter(1)
   MPV_err=res.ParError(1)
   return (MPV, MPV_err)

def Getchi2pNDF(runNr, fixed_ch):
   
   fname=afswork+'Polyparams/run'+str(runNr)+'/polyparams_'+fixed_ch+'.csv'
   if not os.path.exists(fname): return -999.
   with open(fname, 'r') as handle:
      reader=csv.reader(handle)
      alldata=[row for row in reader]
      data=alldata[0]
      chi2pNDF=data.pop() # We don't need to get the chi2 here.
      return float(chi2pNDF)

def Makechi2pNDFdict(runNr, iteration):
   path=afswork+'Polyparams/run'+str(runNr)+'/'
   res={}
   for filename in os.listdir(path):
      fixed_ch=filename[filename.find('2'):filename.find('.csv')]
      
      with open(path+filename, 'r') as handle:
         reader=csv.reader(handle)
         alldata=[row for row in reader]
         if len(alldata)<iteration:
            print(filename)
            continue
         data=alldata[iteration-1]
         chi2pNDF=data.pop()
      res[fixed_ch]=float(chi2pNDF)
   sorted_tuples=sorted(res.items(), key=lambda x:x[1])
   sorted_d={k:v for k,v in sorted_tuples}
   return sorted_d

def GetAttenuationLengthData(runNr, fixed_ch):
   fname=afswork+'attenuationlengths/run'+str(runNr)+'/csvfiles/attenuationlength_'+fixed_ch+'.csv'
   if not os.path.exists(fname): return -999.
   with open(fname, 'r') as handle:
      reader=csv.reader(handle)
      alldata=[row for row in reader]
      data=alldata[0]
      return data

def GetPltThr(func, parameters):
   a, b, c =parameters
   func.SetParameters(a, b, c)
   thr = 0.
   for x in range(0, 6000):
      x_val=x/1000
      y_val = func.Eval(x_val)
      if ROOT.TMath.IsNaN(y_val) or not ROOT.TMath.Finite(y_val):
         thr=x_val
         break
   return x_val

def GetXcalculated(dt, L, cs, wanted=None):
   
   left, right=list(filter(lambda x : x[0]<8, cs)), list(filter(lambda x : x[0]>7, cs))
   NL, NR=len(left), len(right)
   if NL==0 or NR == 0: return -999.
   sumOfInverses=lambda x : sum( [1/i[1] for i in x] )
   A, B = 1/NL*sumOfInverses(left), 1/NR*sumOfInverses(right)
   xcalc = (dt+L*B)/(A+B)

   return xcalc

def Getuncorrectedcscint(runNr, fixed_ch):
   with open(afswork+'uncorrectedcscintvalues/run'+str(runNr)+'/'+fixed_ch+'.csv', r) as h:
      reader=csv.reader(h)
      data=[row for row in reader]
      res=data[0]
   return res

def GetMPV(runNr, fixed_ch, iteration):
   fname=afswork+'MPVs/run'+str(runNr)+'/MPV_'+fixed_ch+'.csv'
   if not os.path.exists(fname): return -999.
   with open(fname, 'r') as h:
      reader=csv.reader(h)
      data=[row for row in reader]
      if data==[]: return -999. # To be investigated
      res=data[iteration]
   MPV, MPV_err = float(res[0]), float(res[1])
   return (MPV, MPV_err)

def MakeFixedCh(fixed):
   fixed_subsystem, fixed_plane, fixed_bar, fixed_SiPM = fixed
   return str(fixed_subsystem)+str(fixed_plane)+'00'+str(fixed_bar)+'_'+str(fixed_SiPM)

def correct_ToF(fixed_SiPM, clock, xpred, cs, xref):
   
   # fixed_subsystem, fixed_plane, fixed_bar, fixed_SiPM = fixed
   time=clock*6.25

   c_SiPM=float(cs[0])
   ToFcorrection=abs((xpred-xref)/c_SiPM)
   if fixed_SiPM<8:
      if xpred >= xref: corrected_t = time - ToFcorrection
      else: corrected_t = time + ToFcorrection
   
   else: 
      if xpred >= xref: corrected_t = time + ToFcorrection
      else: corrected_t = time - ToFcorrection

   return (fixed_SiPM, corrected_t)

def correct_TW(iteration, fixed, rootfile, params, qdc, clock):
      
   highQDCthr=30.
   fixed_ch=MakeFixedCh(fixed)
   plane, bar, SiPM = fixed

   if qdc==-999. or clock==-999.:
      return -999.
   time=clock*6.25

   histname='dtvqdc_'+fixed_ch+'_iteration'+str(iteration)

   # Determine low QDC threshold as highest QDC bin with != 0 entries
   # rootfile e.g. /afs/cern.ch/work/a/aconsnd/Timing/rootfiles/timewalk_dists_24004_4.root
   res=DetermineLowQDCThreshold(rootfile, histname)
   if not isinstance(res, tuple):
      print('Error. result =', res)
      return res
   else:
      a, b =res
   if b==-999.:
      return -999.
   if b==-998.:
      return -998.

   lowQDCthr=b[0]

   if qdc>highQDCthr: corrected_time = time # No correction, append uncorrected times to the final list
   
   elif qdc>lowQDCthr and qdc<=highQDCthr: # Use log parameters computed for iteration i
      A, B, C = params
      TW_correction=A*ROOT.TMath.Log( (B+qdc)/(B+highQDCthr) )
      corrected_time = time + TW_correction
      # corrected_time=(SiPM, tcorr)
   
   else:
      A, B, C = params
      QDCmode, maxEntries=a
      QDC20, max20=b
      dQDC=QDCmode-QDC20
      y2=A*ROOT.TMath.Log(B+QDCmode)+C
      m, c = (y2-max20)/dQDC, y2-(y2-max20)*QDCmode/(dQDC)
      corrected_time = m * qdc + c
      # corrected_time=(SiPM, tcorr) # Investigate efficacy of linear extrapolation method between QDCmax and low threshold

   return (SiPM, corrected_time)

def correct_TW_poly(iteration, fixed, rootfile, polyparams, qdc, clock):
   highQDCthr=30.
   fixed_ch=MakeFixedCh(fixed)
   plane, bar, SiPM = fixed

   if qdc==-999. or clock==-999.:
      return -999.

   histname='dtvqdc_'+fixed_ch+'_iteration'+iteration
   res=DetermineLowQDCThreshold(rootfile, histname)
   if not isinstance(res, tuple):
      return res 

def DetermineLowQDCThreshold(rootfile, histname):
   iteration=rootfile[rootfile.find('iteration')+len('iteration')]

   infile=ROOT.TFile.Open(rootfile, 'read')
   hist=infile.Get(histname) # z.B dtvqdc_24004_4_iteration1

   detID, SiPM = histname.split('_')[1:3]
   lowQDCname='lowQDC_'+detID+'_'+SiPM

   xproj=hist.ProjectionX()
   modalqdc=xproj.GetMaximumBin()
   yproj_modalqdc=hist.ProjectionY('_py', modalqdc, modalqdc)
   meandt_modalqdc=yproj_modalqdc.GetMean()
   maxqdcentries=xproj.GetBinContent(modalqdc)

   a=(modalqdc, meandt_modalqdc)
   binwidth=xproj.GetBinWidth(1)

   b=-998.
   for qdc in range(modalqdc, 0, -1*int(binwidth)):
      tmp=hist.ProjectionY('_tmp'+str(qdc), qdc, qdc)
      tmp_mean=tmp.GetMean()
      if tmp.GetEntries()>0: b=(qdc, tmp_mean)
      else: break

   infile.Close()
   return a, b

def GoodFitFinder(runNr, subsystem, plane, side):
   
   SiPMsdict={'left':(0,1,3,4,5,6), 'right':(8,9,11,12,14,15)}
   chi2s={}
   
   for bar in range(10):
      for SiPM in SiPMsdict[side]:
         fixed_ch = MakeFixedCh((subsystem, plane, bar, SiPM))
         data = GetPolyparams_i(runNr,fixed_ch, 1)
         if data==-999.:continue
         chi2pNDF=data[-1]
         chi2s[chi2pNDF]=(bar, SiPM)
            
   return chi2s

def GetCutDistributions(runNr, distmodes):
   Allmodes=('yresidual', 'nSiPMs', 'slopes')
   datafile=afswork+'SelectionCriteria/SelectionCriteria_run'+runNr+'.root'

   if isinstance(distmodes, str):
      distmodes=(distmodes,)

   for distmode in distmodes:
      if distmode not in Allmodes:
         print('Use a valid mode.')
         return 0

   file=ROOT.TFile.Open(datafile, 'READ')
   dists={}

   for distmode in distmodes:
      if distmode=='yresidual':   
         for s in (1,2):
            for p in range(systemAndPlanes[s]):
               key=str(s*10+p)
               name=key+'_yresidual'
               hist=file.Get(name).Clone()
               hist.SetDirectory(ROOT.gROOT)
               dists[name]=hist

      else:
         hist=file.Get(distmode).Clone()
         hist.SetDirectory(ROOT.gROOT)
         dists[distmode]=hist

   return dists
