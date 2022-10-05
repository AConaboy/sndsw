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

def GetSide(fixed_ch):
    detID, SiPM = fixed_ch.split('_')
    s,p,b=parseDetID(int(detID))
    if s==2:
        if int(SiPM)<8: side='left'
        elif int(SiPM)>7: side='right'
    elif s==1:
        if int(SiPM)<6: side='left'
        elif int(SiPM)>5: side='right'
    return side

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

def OneHitPerSystem(hits, systems):
    #USPlanes={k:0 for k in range(5)}
    #hitdict={int(k)+10*int(s):0 for k in range(systemAndPlanes[s]) for s in systems}
    hitdict={}
    for s in systems:
        for k in range(systemAndPlanes[s]):
            key=10*s+k
            hitdict[key]=0
    for i, hit in enumerate(hits):
        if not hit.isValid(): continue
        detID=hit.GetDetectorID()
        s, p, b = parseDetID(detID)
        if s==3: continue
        key=10*s+p
        hitdict[key]+=1

    for key in hitdict:
        if hitdict[key] != 1: return False
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
    if s==1: 
        nSiPMs, SiPMs_plane=12, 84 # is it? 
        return SiPM+nSiPMs*b+p*SiPMs_plane
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
    if s == 1: return SiPM*b*12
    if s == 2:  return SiPM+b*16 

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

def GetDSH_average(hits, nPlanes=3):
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

def Getcscint(runNr, fixed_ch, iteration):
    if not os.path.exists(afswork+'cscintvalues/run'+str(runNr)+'/cscint_'+fixed_ch+'.csv'): return -999.
    with open(afswork+'cscintvalues/run'+str(runNr)+'/cscint_'+fixed_ch+'.csv', 'r') as handle:
        reader=csv.reader(handle)
        alldata=[row for row in reader]
        if len(alldata)<int(iteration) or len(alldata) == 0: return -999.
        else:
            data=alldata[int(iteration)]
            return (float(data[1]), float(data[2]))

def Getcscint_chi2pNDF(runNr,fixed_ch,iteration):
   # with open(afswork+'rootfiles/run'+str(runNr)+'/cscintvalues/cscint_'+fixed_ch+'.csv', 'r') as handle:
   with open(afswork+'cscintvalues/run'+str(runNr)+'/cscint_'+fixed_ch+'.csv', 'r') as handle:
      reader=csv.reader(handle)
      alldata=[row for row in reader]
      if len(alldata)<iteration or len(alldata) == 0: return -999.
      else:
         data=alldata[iteration]
         return float(data[-2])/int(data[-1])

def Getcscintdict(runNr, subsystem, iteration):
   path=f'{afswork}/cscintvalues/run{runNr}/'
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
         data=alldata[iteration]
      res[fixed_ch]=float(data[1])
   sorted_tuples=sorted(res.items(), key=lambda x:x[1])
   sorted_d={k:v for k,v in sorted_tuples}
   return sorted_d
   
def GetAvgcscint(path, detID, SiPMs, iteration):
   cs=[]
   for SiPM in SiPMs:
      fixed_ch=('_').join( (str(detID), str(SiPM)) ) 
      val=Getcscint_i(path, fixed_ch, iteration)
      cs.append(val)
   
   cbar = sum([i[0] for i in cs])/len(cs)
   cbar_err = ROOT.TMath.Sqrt( sum( [i[1]**2 for i in cs] ) )
   return cbar, cbar_err

def GetPolyParams(runNr, fixed_ch, n, iteration):
    
   filelengths={1:11,2:13,3:15,5:9}
   if not os.path.exists(f'{afswork}Polyparams/run{runNr}/polyparams{n}_{fixed_ch}.csv'): return -999.
   with open(f'{afswork}Polyparams/run{runNr}/polyparams{n}_{fixed_ch}.csv', 'r') as f:
      reader=csv.reader(f)
      alldata=[r for r in reader]
      if alldata==[]:return -999
      #print(f'All data {fixed_ch}, len:{len(alldata[iteration-1])}: {alldata}')
      data=alldata[int(iteration) - 1]
   # print('=='*15,f'\nlen(data):{len(data)}\n','=='*15)
   if len(data)!=filelengths[n]: return -999
     
   params=[float(i) for i in data[1:-2]]
   limits=[float(i) for i in data[-2:]]
   return params,limits

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
   fname=f'{afswork}rootfiles/run{runNr}/timewalk_{fixed_ch}.root'
   if not os.path.exists(fname): return -999.
   f=ROOT.TFile.Open(fname, 'READ')
   histname=f'dtvqdc_{fixed_ch}_iteration{iteration}'
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

def Getchi2_info(runNr, fixed_ch, n, iteration):
    fname=f'{afswork}chi2s/run{runNr}/chi2s{n}_{fixed_ch}.csv'
    if not os.path.exists(fname): return -999.
    with open(fname, 'r') as handle:
        reader=csv.reader(handle)
        alldata=[row for row in reader]
        if len(alldata)==0 or len(alldata)<iteration: return -999.
        data=alldata[iteration-1]
        try:
            chi2info=int(data[0]), float(data[1]), int(data[2])
        except ValueError:
            print(f'Non-integer NDF for {fixed_ch}')
            return -999.
    return chi2info

def Getchi2pNDF(runNr, fixed_ch, n, iteration):
    fname=f'{afswork}chi2s/run{runNr}/chi2s{n}_{fixed_ch}.csv'
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

def GetNDF(runNr, fixed_ch, iteration):
    fname=afswork+'chi2s/run'+str(runNr)+'/chi2s_'+fixed_ch+'.csv'
    if not os.path.exists(fname): return -999.
    with open(fname, 'r') as handle:
        reader=csv.reader(handle)
        alldata=[row for row in reader]
        if len(alldata)==0 or len(alldata)<iteration: return -999.
        data=alldata[iteration-1]
    return int(data[2])

def GetBadchi2pNDFdict(runNr, subsystem, iteration):
   #path=afswork+'chi2s/run'+str(runNr)+'/'
   path=f'{afswork}/chi2s/run{runNr}/'
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

def Makechi2pNDFdict(runNr, subsystem, iteration):
   #path=afswork+'chi2s/run'+str(runNr)+'/'
   path=f'{afswork}/chi2s/run{runNr}/'
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
   fname=f'{afswork}MPVs/run{runNr}/MPV_{fixed_ch}.csv'
   if not os.path.exists(fname): return -999.
   with open(fname, 'r') as h:
      reader=csv.reader(h)
      data=[row for row in reader]
      if data==[]: return -999. # To be investigated
      res=data[iteration-1]
   return float(res[0])

def MakeFixedCh(fixed):
   fixed_subsystem, fixed_plane, fixed_bar, fixed_SiPM = fixed
   return str(fixed_subsystem)+str(fixed_plane)+'00'+str(fixed_bar)+'_'+str(fixed_SiPM)

def correct_ToF(SiPM, clock, xpred, cs, xref):
   
   # fixed_subsystem, fixed_plane, fixed_bar, fixed_SiPM = fixed
   time=clock*6.25

   c_SiPM=float(cs[0])
   ToFcorrection=abs((xpred-xref)/c_SiPM)
   if SiPM<8:
      if xpred >= xref: corrected_t = time - ToFcorrection
      else: corrected_t = time + ToFcorrection
   
   else: 
      if xpred >= xref: corrected_t = time + ToFcorrection
      else: corrected_t = time - ToFcorrection

   return (SiPM, corrected_t)

def GoodFitFinder(runNr, subsystem, plane, side):
   
   SiPMsdict={'left':(0,1,3,4,6,7), 'right':(8,9,11,12,14,15)}
   chi2s={}
   
   for bar in range(10):
      for SiPM in SiPMsdict[side]:
         fixed_ch = MakeFixedCh((subsystem, plane, bar, SiPM))
         chi2pNDF = Getchi2pNDF(runNr,fixed_ch, 1)
         if chi2pNDF==-999.:continue
         #chi2pNDF=data[-1]
         chi2s[chi2pNDF]=(bar, SiPM)
            
   return chi2s

def GetCutDistributions(runNr, distmodes):
   Allmodes=('yresidual', 'nSiPMs', 'slopes')
   filename=f'{afswork}SelectionCriteria/SelectionCriteria_run{runNr}.root'
   if not os.path.exists(filename): filename=f'{afswork}SelectionCriteria/SelectionCriteria_run004612.root'

   if isinstance(distmodes, str):
      distmodes=(distmodes,)

   for distmode in distmodes:
      if distmode not in Allmodes:
         print('Use a valid mode.')
         return 0

   f=ROOT.TFile.Open(filename, 'READ')
   dists={}

   for distmode in distmodes:
      if distmode=='yresidual':   
         for s in (1,2):
            for p in range(systemAndPlanes[s]):
               key=str(s*10+p)
               name=key+'_yresidual'
               hist=f.Get(name).Clone()
               hist.SetDirectory(ROOT.gROOT)
               dists[name]=hist

      else:
         hist=f.Get(distmode).Clone()
         hist.SetDirectory(ROOT.gROOT)
         dists[distmode]=hist

   return dists

def adjacentplaneschi2s(runNr, subsystem, plane, side):
    if plane != 4 and plane != 0: to_try=(plane-1, plane+1)
    elif plane==0: to_try=(1,2)
    elif plane==4: to_try=(2,3)

    best_efforts={}

    for testplane in to_try:
        chi2s=GoodFitFinder(runNr, subsystem, plane, side)
        if len(chi2s.keys())==0: continue
        tmp=min(list(chi2s.keys()))
        best_efforts[tmp]=plane

    return best_efforts

def GetEntriesInHist(runNr, fixed_ch, mode, iteration):
    f=ROOT.TFile.Open(f'{afswork}/rootfiles/run{runNr}/timewalk_{fixed_ch}.root', 'READ')
    hist=f.Get(f'{mode}_{fixed_ch}_iteration{iteration}')
    hist.SetDirectory(ROOT.gROOT)
    entries=hist.GetEntries()
    f.Close()
    return entries

def GetPolyParamRanges(runNr, n, subsystem, iteration):
    params={i:[0.,0.] for i in ('A', 'B', 'C')}
    path=f'{afswork}Polyparams/run{runNr}/'
    files=[i for i in os.listdir(path) if int(i.split('_')[1][0])==subsystem]
    counter=0
    for f in files:
        if f.split('_')[0].find(str(n)) == -1:continue
        #if idx>1: break
        fixed_ch=f"{f.split('_')[1]}_{f.split('.')[0].split('_')[-1]}"
        tmp = GetPolyParams(runNr, fixed_ch, n, iteration)
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

def MakePolyParamDicts(runNr, subsystem, iteration):
    #A,B,C={}, {}, {}
    params={i:[] for i in ('A', 'B', 'C')}
    path=f'{afswork}Polyparams/run{runNr}/'
    files=[i for i in os.listdir(path) if int(i.split('_')[1][0])==subsystem]
    for idx, f in enumerate(files):
        fixed_ch=f"{f.split('_')[1]}_{f.split('.')[0].split('_')[-1]}"
        fps=GetPolyParams(runNr, fixed_ch, iteration)
        if fps==-999.: continue
        vals=[fps[i*2] for i in range(3)]
        [params[p].append(vals[i]) for i,p in enumerate(params)]
    ranges={i:[0,0] for i in ('A', 'B', 'C')}
    #[ranges[i]=[min(params[i]), max(params[i])] for i in ('A', 'B', 'C')]
    #return params, ranges
    return params

def checkFitStatus(fitter):
    for attribute in ('hist', 'graph'): 
        if isinstance(getattr(fitter, attribute), float)==True: return False
    return True

   
