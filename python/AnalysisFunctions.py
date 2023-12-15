#!/usr/bin/env python
import ROOT,os,csv,json
import math as m 

# class mufi_analysis: # I could make several classes for generic functions like parseDetID and more specific ones for 

class Analysis(object):

	def __init__(self, options):
		self.options=options

		# Adding flag for LaserMeasurements/ work to use some functions defined in here
		if not hasattr(options, "LaserMeasurements"):
			if not hasattr(options, 'runNumber'): options.runNumber=options.runs[0]
			self.runNr = str(options.runNumber).zfill(6)
			self.TWCorrectionRun = str(options.TWCorrectionRun).zfill(6)
			self.freq=160.316E6
			self.TDC2ns=1E9/self.freq
			self.timealignment=self.GetTimeAlignmentType(self.runNr)
			self.state=options.state
			if hasattr(options, 'datafiletype'): self.fileext=options.datafiletype
			else: self.fileext='csv'		
			self.CorrectionType=options.CorrectionType

			afswork='/afs/cern.ch/work/a/aconsnd/Timing'
			afsuser='/afs/cern.ch/user/a/aconsnd/twfiles'
			# if self.options.path.find('commissioning/TI18')>0:
			if options.datalocation=='commissioning':
				self.path=afswork+'-commissioning/'
			elif options.datalocation=='physics':
				self.path=afswork+'-physics2022/'
			elif options.datalocation=='H8':
				self.path=afswork+'-H8/'

		self.correctionparams=lambda ps : [y for x,y in enumerate(ps) if x%2==0]
		self.correctionfunction = lambda ps, qdc : ps[3]*(qdc-ps[0])/( ps[1] + ps[2]*(qdc-ps[0])*(qdc-ps[0]) ) + ps[4]*(qdc-ps[0])
		self.A, self.B = ROOT.TVector3(), ROOT.TVector3()
		self.systemAndPlanes = {1:2,2:5,3:7}
		self.systemAndBars = {1:7,2:10,3:60}
		self.systemAndChannels = {1:[0,8],2:[2,6],3:[0,1]}
		self.systemAndSiPMs={1:range(16),2:(0,1,3,4,6,7,8,9,11,12,14,15),3:(1,)}
		self.verticalBarDict={0:1, 1:3, 2:5, 3:6}
		self.gelsides={0:'right', 1:'left', 2:'right', 3:'left', 4:'left'}
		self.subsystemNames={1:'veto', 2:'upstream', 3:'downstream'}
		self.verbose=False
		
		if hasattr(options, 'datafiletype'): self.fileext=options.datafiletype
		else: self.fileext='csv'

		self.sigmatds0=0.263, 9.5E-5

	def SetTask(self, task):
		self.task=task

	def DSHcheck(self, detID): # True if detID is a horizontal bar
		s,p,b=self.parseDetID(detID)
		if s!=3: return False 
		if p<3 and b<60:  return True
		else: return False

	def DSVcheck(self, detID): # True if detID is a vertical 
		s,p,b=self.parseDetID(detID)
		if s!=3: return False
		if p<3 and b>59: return True 
		elif p==3: return True 
		else: return False

	def GetListOfChannels(self, subsystem): #Only returns horizontal channels for the DS (s==3)
		channels=[f'{self.MakeFixedCh((subsystem, plane, bar, SiPM))}' for plane in range(self.systemAndPlanes[subsystem]) for bar in range(self.systemAndBars[subsystem]) for SiPM in self.systemAndSiPMs[subsystem] ]
		return channels

	def GetGeoFile(self,runNumber):
		if type(runNumber)==str: runNumber=int(runNumber)
		
		if runNumber < 4575: geofile='V3_08August2022'
		elif runNumber < 4856: geofile='V5_14August2022'
		elif runNumber < 4992: geofile='V6_08October2022'
		elif runNumber > 4991: geofile='V7_22November2022'

		return f"geofile_sndlhc_TI18_{geofile}.root"	

	def GetRunYear(self, runNr):
		if isinstance(runNr, str): runNr=int(runNr)

		if runNr < 5485: year='2022'
		else: year='2023'

		return year

	def BuildBarLengths(self, MuFilter):
		Vetobarlength = MuFilter.GetConfParF('MuFilter/VetoBarX')
		USbarlength = MuFilter.GetConfParF('MuFilter/UpstreamBarX')
		DSbarlength_hor = MuFilter.GetConfParF('MuFilter/DownstreamBarX')
		DSbarlength_vert = MuFilter.GetConfParF('MuFilter/UpstreamBarY_ver')
		barlengths={1:Vetobarlength, 2:USbarlength, 3:DSbarlength_hor, 4:DSbarlength_vert}
		return barlengths

	def BuildzPos(self, MuFilter, Scifi):
		A,B=ROOT.TVector3(), ROOT.TVector3()
		zPos={'MuFilter':{},'Scifi':{}}
		# MuFilter
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
		# Scifi
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
		if s<3:
			if int(SiPM)<8: side='left'
			elif int(SiPM)>7: side='right'
			else: return
		elif s==3:
			s,p,b=self.parseDetID(int(detID))
			if p!=3 and b<60 and int(SiPM)==0: side='left'
			elif p!=3 and b<60 and int(SiPM)==1: side='right'
			elif p!=3 and b>59 and int(SiPM)==0: side='top'
			elif p==3 and int(SiPM)==0: side='top'
			else: 
				print('huh?')
				return

		return side

	def IsGel(self, fixed_ch):
		detID, SiPM = fixed_ch.split('_')
		s,p,b = self.parseDetID(int(detID))
		if s!=2: return 0
		if int(SiPM)>7: side='right'
		else: side='left'
		if self.gelsides[p]==side: return 1
		else: return 0
  
	def parseDetID(self, detID):
		if not isinstance(detID, int): detID=int(detID)
		subsystem=detID//10000
		plane=detID%10000//1000
		bar=detID%1000
		return subsystem, plane, bar

	def MakeDetID(self, fixed):
		fixed_subsystem, fixed_plane, fixed_bar = fixed
		if fixed_subsystem in (1,2,3):
			return int(f'{str(fixed_subsystem)}{str(fixed_plane)}{str(fixed_bar).zfill(3)}')

	def MakeFixedCh(self, fixed):
		fixed_subsystem, fixed_plane, fixed_bar, fixed_SiPM = fixed
		return f'{str(fixed_subsystem)}{str(fixed_plane)}{str(fixed_bar).zfill(3)}_{str(fixed_SiPM)}'

	def GetDSPlaneNumber(self, detID):
		s,p,b=self.parseDetID(detID)
		if s<3: return p 
		if b>59 and p<3: plane=p*2+1
		elif any( [b<60, b>59 and p==3] ): plane=p*2 
		return plane

	def MakeHumanReadableFixedCh(self, fixed_ch):
		s,p,b=self.parseDetID(int(fixed_ch.split('_')[0])) 
		SiPM=int(fixed_ch.split('_')[1])
  		# Not using f-string because ROOT can't cope with combining them and TLatex
		res=str(self.subsystemNames[s])+', plane '+str(p+1)+', bar '+str(b+1)+', SiPM '+str(SiPM+1)+', '+self.GetSide(fixed_ch)+' side' # +1 to make things more readable. I hope this doesn't complicate things
		return res

	def MakeHumanReadableDetID(self, detID):
		s,p,b=self.parseDetID(detID) 
  		# Not using f-string because ROOT can't cope with combining them and TLatex
		res=str(self.subsystemNames[s])+', plane '+str(p+1)+', bar '+str(b+1)
		return res

	def GetScifiAverageTime(self, hits):
		pass

	def GetDSHaverage(self, hits, mode='tds0'):
		stations={k:{} for k in range(4)} # only k=0,1,2 are used
		for i,hit in enumerate(hits):
			detID=hit.GetDetectorID()
			s,p,b=self.parseDetID(detID)
			if s!=3:continue
			if b>59: continue # The bar IDs in DS4V start at 60 so the vertical bars are not contaminating here.
			
			stations[p][i]=hit

		# if not all( [len(stations[p])==1 for p in range(nPlanes)] ): return -999. # Event selection: only accept events with 1 bar in all DS planes
		total=0.
		counter=0
		if mode=='tds0':
			for p in range(4):
				hits=list(stations[p].values())
				for hit in hits:
					if not hit.isValid():continue
					detID=hit.GetDetectorID()
					tdcs=hit.GetAllTimes()
					if len(tdcs) != 2: continue # pass
					# for all 3 horizontal DS planes.
					for item in tdcs: 
						SiPM, clock = item
						dscorrectedtime=self.task.MuFilter.GetCorrectedTime(detID, SiPM, clock*self.TDC2ns, 0)
						total+=dscorrectedtime
						counter+=1
			if counter==0: 
				return -999

			dsh_average=total/counter
			
			return dsh_average, counter

		elif mode=='timingdiscriminant': # Take k=2 for the 3rd horizontal plane
			DS3H_hits=list(stations[2].values())
			total=[]

			for DS3H_hit in DS3H_hits:
				if not DS3H_hit.isValid(): continue
				detID=DS3H_hit.GetDetectorID()
				times=DS3H_hit.GetAllTimes()
				if len(times)!=2: continue # Require left and right SiPMs fire 
				
				for idx in times: 
					SiPM,clock=idx
					dscorrectedtime=self.task.MuFilter.GetCorrectedTime(detID, SiPM, clock*self.TDC2ns, 0)
					total.append(dscorrectedtime)
			if len(total)==0:
				return -420					
			return sum(total)/len(total)

		elif mode=='deltastations':
			res={f'delta32':0, 'delta21':0}
			ts={k:[] for k in range(3)}
			for i in range(3):
				dshits=list(stations[i].values())
				
				for hit in dshits:
					detID=hit.GetDetectorID()
					times=hit.GetAllTimes()

					for idx in times:
						SiPM, clock=idx
						dscorrectedtime=self.task.MuFilter.GetCorrectedTime(detID, SiPM, clock*self.TDC2ns, 0)
						ts[i].append(dscorrectedtime)
			if any([len(ts[i])==0 for i in ts]): return False
			res['delta32']=sum(ts[2])/len(ts[2]) - sum(ts[1])/len(ts[1])
			res['delta21']=sum(ts[2])/len(ts[2]) - sum(ts[0])/len(ts[0])
			return res

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

	def GetAverageTime(self, mufiHit, side='both'):
		value=[0, 0]
		count=[0, 0]
		nSiPMs=mufiHit.GetnSiPMs()
		times=mufiHit.GetAllTimes()
		
		s, p, b = self.parseDetID(mufiHit.GetDetectorID())
		for element in times:
			SiPM, clock = element
			if s==2 and self.IsSmallSiPMchannel(SiPM):continue
			if SiPM<nSiPMs:
				value[0]+=clock*self.TDC2ns
				count[0]+=1
			else:
				value[1]+=clock*self.TDC2ns
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
			else: 
				return
				# print(f'hit does not have fired SiPMs on both sides.')
		elif s == 3 and b < 60: # I don't require a SiPM to fire on both sides here because that is required by the TDS0 determination
			if side=='both':
				average = 0.5*(value[0]/count[0]+value[1]/count[1])
				return average
			elif side == 'L': return value[0]/count[0]
			elif side == 'R': return value[1]/count[1]
   
	# def GetMedianTime(self, hit, mode=None, particleToF=None):
	def GetMedianTime(self, hit, mode=None):
		
		if not mode: mode=self.state
		# if self.timealignment=='old' and particleToF==None: 
		# 	print(f'A particle ToF must be provided for old time alignment data')
		# 	return -999

		values={'left':[], 'right':[]}
		nSiPMs=hit.GetnSiPMs()
		clocks=hit.GetAllTimes()
		if mode != 'uncorrected': qdcs=hit.GetAllSignals()
		detID=hit.GetDetectorID()
		s, p, b = self.parseDetID(detID)
		
		# Should not be needed
		if s==3 and b>59: 
			vals=[i[1] for i in clocks]
			return sum(vals)/len(vals)
  
		nLeft, nRight = self.GetnFiredSiPMs(hit)
		nSiPMconditions={1: any( [nLeft<6, nRight<6]),
						2: any( [nLeft<4, nRight<4]),
						3: all(( b<60, any( [nLeft!=1, nRight!=1]) ))
      					}
		if nSiPMconditions[s]: 
			# print(f'nSiPM conditions not met for {detID}')
			return
			
		for element in clocks:
			SiPM, clock = element
			fixed_ch=self.MakeFixedCh((s,p,b,SiPM))		
			side=self.GetSide(fixed_ch)

			if s==2 and self.IsSmallSiPMchannel(SiPM): continue

			if mode == 'uncorrected': values[side].append(clock*self.TDC2ns)

			elif mode == 'corrected':
				qdc=self.GetChannelVal(SiPM, qdcs)
				if not qdc: continue
				correctedtime=self.correct_TW(fixed_ch, qdc, clock)
				if not correctedtime: continue
				values[side].append(correctedtime)
			
			elif mode == 'alignment':
				qdc=self.GetChannelVal(SiPM, qdcs)
				if not qdc: continue
				correctedtime=self.MuFilterCorrectedTime(fixed_ch, qdc, clock)
				if not correctedtime: continue
				
				values[side].append(correctedtime)
	
		if not all( (len(values[i])>0 for i in values) ): 
			print(f'Zero entries on one or both sides of bar {detID}')
			return
		medians={}
  
		if s==3 and b<60 and all( [nLeft==1, nRight==1] ): 
			for x in values: medians[x]=values[x][0]
			return medians

		for x in values:
			if len(values[x])%2==0:
				if len(values[x])==2: print(detID, values)
				medians[x] = 0.5* ( values[x][int(len(values[x])/2)] + values[x][int( len(values[x])/2+1)] )
			else: 
				medians[x] = values[x][int( 0.5*(len(values[x])+1))]
		return medians


	def GetChannelVal(self, SiPM, chs):
		for entry in chs:
			fSiPM, val = entry
			if fSiPM == SiPM:
				return val
		return

	def OneHitPerSystem(self, hits, systems, Nfired=False):
		verbose=self.verbose

		hitdict={}
		for s in systems: # systems always includes US and DS. Veto is included if a Scifi track is also formed.
			for p in range(self.systemAndPlanes[s]):
				key=10*s+p
				hitdict[key]=[]
		for hit in hits:
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
			if key<30 and len(hitdict[key]) != 1: return False

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

	def GetSubsystemZlimits(self, subsystem, zPos):
		if subsystem==0:
			pass
		elif subsystem==1:
			zmin, zmax=self.zPos['MuFilter'][10], self.zPos['MuFilter'][11]
		elif subsystem==2:
			zmin, zmax=self.zPos['MuFilter'][20], self.zPos['MuFilter'][24]
		elif subsystem==3:
			zmin, zmax=self.zPos['MuFilter'][30], self.zPos['MuFilter'][36]
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

	def GetnFiredSiPMs(self, hit):
		
		s,p,b=self.parseDetID(hit.GetDetectorID())
		nSiPMs=hit.GetnSiPMs()
		
		nFiredSiPMs_left=0
		nFiredSiPMs_right=0
		channels=hit.GetAllSignals()
		for ch in channels:
			SiPM, qdc = ch 
			if s==2 and self.IsSmallSiPMchannel(SiPM):continue
			if SiPM<nSiPMs: nFiredSiPMs_left+=1
			elif SiPM>=nSiPMs: nFiredSiPMs_right+=1
		return nFiredSiPMs_left, nFiredSiPMs_right 

	def GetFiredSiPMsOnPCBs(self, hits):
		channels_on_PCB = {}
		for hit in hits:
			detID = hit.GetDetectorID()
			s,p,b = self.parseDetID(detID)
			Fired_left, Fired_right = self.GetnFiredSiPMs(hit)
			fired = {'left':Fired_left, 'right':Fired_right}

			for side in fired: 
				key=f'{p}_{side}'
				if not key in channels_on_PCB: channels_on_PCB[key]=0 
				channels_on_PCB[key] += fired[side]
   
		return channels_on_PCB

	def GetSiPMNumberInSystem_LandR(self, detID, SiPM): # 20000 SiPM 8 -> 8
		if not isinstance(SiPM, int): SiPM=int(SiPM)
		s, p, b = self.parseDetID(int(detID))
		if s==1:
			nSiPMs, SiPMs_plane=16, 112 
			return int(SiPM)+nSiPMs*b+p*SiPMs_plane
		elif s==2:
			nSiPMs, SiPMs_plane=16, 160
			return SiPM+nSiPMs*b+p*SiPMs_plane

		elif s==3: # Count left and right horizontal SiPMs consecutively
			nSiPMs_bar_hor, nSiPMs_bar_ver=2, 1
			nSiPMs_plane_hor, nSiPMs_plane_ver=120, 60
			if b>59 and p<3: # First 3 vertical layers
				horizontalSiPMs=(p+1)*nSiPMs_plane_hor
				verticalSiPMs=p*nSiPMs_plane_ver
				total=horizontalSiPMs+verticalSiPMs+(b-60)*nSiPMs_bar_ver
			elif b<60: # All horizontal layers
				horizontalSiPMs=p*nSiPMs_plane_hor
				verticalSiPMs=p*nSiPMs_plane_ver
				total=horizontalSiPMs+verticalSiPMs+b*nSiPMs_bar_hor+SiPM
			elif b>59 and p==3: # Final vertical layer
				total=3*120+3*60+(bar-60)
		return total

	def GetSiPMNumberInSystem_PCBbyPCB(self, detID, SiPM):
		if not isinstance(detID, int): detID=int(detID)
		if not isinstance(SiPM, int): SiPM=int(SiPM)

		s,p,b=self.parseDetID(detID)
		if s==2:
			nSiPMs=8
			SiPMs_plane=160
			total=b*nSiPMs+p*SiPMs_plane
			if SiPM<nSiPMs:
				return int(total+SiPM )
			else:
				return int(total+SiPMs_plane/2+(SiPM-nSiPMs))
		elif s==1:
			nSiPMs=8
			SiPMs_plane=112
			total=b*nSiPMs+p*SiPMs_plane
			if SiPM<nSiPMs:
				return int(total+SiPM )
			else:
				return int(total+SiPMs_plane/2+(SiPM-nSiPMs))

		elif s==3:
			nSiPMs_bar_hor, nSiPMs_bar_ver=2,1
			nSiPMs_plane_hor, nSiPMs_plane_ver=120, 60
			if b>59 and p!=3: # First 3 vertical layers
				horizontalSiPMs=int((p+1)*nSiPMs_plane_hor)
				verticalSiPMs=int(p*nSiPMs_plane_ver)
				total=horizontalSiPMs+verticalSiPMs+int((b-60)*nSiPMs_bar_ver)
				return total
			elif b<60: # All horizontal layers
				horizontalSiPMs=int(p*nSiPMs_plane_hor)
				verticalSiPMs=int(p*nSiPMs_plane_ver)
				if SiPM==0:
					total=horizontalSiPMs+verticalSiPMs+b
					return total 
				elif SiPM==1:
					total=horizontalSiPMs+verticalSiPMs+b+60
					return total
			elif b>59 and p==3:
				total=3*120+3*60+(b-60)
				return int(total)

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

	def SiPM2BarAndPosition(self, SiPM):
		# Pass the SiPM number on the PCB 
		# returns the bar number and SiPM number within the bar
		barNumber = (SiPM-1)//8+1
		SiPMNumber = (SiPM-1)%8+1
		if SiPMNumber==0: SiPMNumber==8
		return barNumber, SiPMNumber

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
		if count[0] != 0 and count[1] != 0:
			return (mean[0]/count[0]-mean[1]/count[1])/2.
		else: return	

	def Getyresidual(self, detID):
		s,p,b=self.parseDetID(detID)
		key=10*s+p
		z=self.zPos['MuFilter'][key]
		lam=(z-self.task.pos.z())/self.task.mom.z()
		pq = self.task.A-self.task.pos
		uCrossv= (self.task.B-self.task.A).Cross(self.task.mom)
		doca = pq.Dot(uCrossv)/uCrossv.Mag()
		return doca

	def GetdtCalc(self, xpred, L, cs):
		left, right=list(filter(lambda x : x[0]<8, cs)), list(filter(lambda x : x[0]>7, cs))
		NL, NR=len(left), len(right)
		if NL == 0 or NR == 0: return 
		sumOfInverses=lambda x : sum( [1/i[1] for i in x] )
		return xpred/NL*sumOfInverses(left) - (L-xpred)/NR*sumOfInverses(right)

	def Getcscint(self, runNr, fixed_ch, state):

		mode=self.fileext

		if mode=='csv':
			iteration=0 if state=='uncorrected' else 1
			if not os.path.exists(f'{self.path}cscintvalues/run{runNr}/cscint_{fixed_ch}.csv'): 
				return
			with open(f'{self.path}cscintvalues/run{runNr}/cscint_{fixed_ch}.csv', 'r') as handle:
				reader=csv.reader(handle)
				alldata=[row for row in reader]
				
				if len(alldata)<iteration+1 or len(alldata) == 0: return 
				try: data=alldata[iteration]
				except IndexError: print(f'{fixed_ch} IndexError')

			return (float(data[1]), float(data[2]))

		elif mode=='json':
			filename=f'{self.path}cscintvalues/run{runNr}/cscint_{fixed_ch}.json'

			if not os.path.exists(filename): return 

			with open(filename, 'r') as x:
				d=json.load(x)

			if state not in d: return 

			return d[state][0], d[state][1]

	def Getcscint_offset(self, runNr, fixed_ch, state):

		filename=f'{self.path}cscintvalues/run{runNr}/cscint_{fixed_ch}.json'

		if not os.path.exists(filename): return 

		with open(filename, 'r') as x:
			d=json.load(x)

		if state not in d: return 

		return d[state][2], d[state][3]		


	def GetBarAveragecscint(self, runNr, detID, state):

		subsystem, plane, bar = self.parseDetID(detID)
		cscintvalues={}
		for SiPM in self.systemAndSiPMs[subsystem]:
			fixed_ch=self.MakeFixedCh((subsystem, plane, bar, SiPM))
			if fixed_ch not in self.cscintvalues:continue
			cscint=self.cscintvalues[fixed_ch]
			if not cscint: continue
			cscintvalues[fixed_ch]=cscint

		average_cscint = sum( [cscintvalues[ch][0] for ch in cscintvalues] ) / len(cscintvalues.items())
		uncertainty_sq = sum( [cscintvalues[ch][1]**2 for ch in cscintvalues] ) 
		uncertainty=ROOT.TMath.Sqrt(uncertainty_sq)
	
		return average_cscint, uncertainty

	def GetBarAveragesigmat(self, runNr, detID, state):

		subsystem, plane, bar = self.parseDetID(detID)
		bartimeresolutionvalues={}
		for SiPM in self.systemAndSiPMs[subsystem]:
			fixed_ch=self.MakeFixedCh((subsystem, plane, bar, SiPM))
			if fixed_ch not in self.timeresolutionvalues:continue
			timeresolutionvalue=self.timeresolutionvalues[fixed_ch]
			if m.isnan(timeresolutionvalue[0]): continue

			bartimeresolutionvalues[fixed_ch]=timeresolutionvalue

		average_timeresolution = sum( [bartimeresolutionvalues[ch][0] for ch in bartimeresolutionvalues] ) / len(bartimeresolutionvalues.items())
		uncertainty_sq = sum( [bartimeresolutionvalues[ch][1]**2 for ch in bartimeresolutionvalues] ) 
		uncertainty=ROOT.TMath.Sqrt(uncertainty_sq)
	
		return average_timeresolution, uncertainty	

	# def GetBarAverageQDC(self, qdcs):
		# for 

	def Getcscint_chi2pNDF_info(self, runNr,fixed_ch,state):
		iteration=0 if state=='uncorrected' else 1
		with open(f'{self.path}cscintvalues/run{runNr}/cscint_{fixed_ch}.csv', 'r') as handle:
			reader=csv.reader(handle)
			alldata=[row for row in reader]
			if len(alldata)<iteration+1: return 
			try:
				data=alldata[iteration]
			except IndexError:
				print(f'{fixed_ch} IndexError')
		return float(data[-2]),int(data[-1])

	def Getcscint_chi2pNDF(self, runNr,fixed_ch,state):
		iteration=0 if state=='uncorrected' else 1
		if not os.path.exists(f'{self.path}cscintvalues/run{runNr}/cscint_{fixed_ch}.csv'): return
		with open(f'{self.path}cscintvalues/run{runNr}/cscint_{fixed_ch}.csv', 'r') as handle:
			reader=csv.reader(handle)
			alldata=[row for row in reader]
			if len(alldata)<iteration+1: return 
			try:
				data=alldata[iteration]
			except IndexError:
				print(f'{fixed_ch} IndexError')
		return float(data[-2])/int(data[-1])

	def Makecscintdict(self, runNr, state):
		d={}
		for s in (1,2,3):
			for p in range(self.systemAndPlanes[s]):
				for b in range(self.systemAndBars[s]):
					for SiPM in self.systemAndSiPMs[s]:
						fixed_ch=self.MakeFixedCh((s,p,b,SiPM))
						cscint=self.Getcscint(runNr, fixed_ch=fixed_ch, state=state)
						if not cscint: continue
						else: d[fixed_ch]=cscint
		if len(d)==0: self.cscintvalues=None
		self.cscintvalues=d

	def Maketimeresolutiondict(self, runNr, state):
		sigmadelta=2/ROOT.TMath.Sqrt(6)*self.sigmatds0[0]
		iteration=0 if state=='uncorrected' else 1
		path=f'{self.path}TimeResolution/run{runNr}/'
		res={}
		for filename in os.listdir(path):
			fixed_ch=f'{filename.split("_")[1]}_{filename.split("_")[2].split(".")[0]}'
			if fixed_ch[0]==str(3): continue
			with open(path+filename, 'r') as handle:
				reader=csv.reader(handle)
				alldata=[row for row in reader]
				if len(alldata)<iteration or alldata==[]:
					continue
				data=alldata[iteration]
				# data -> (state, time resolution, uncertainty on time resolution) #
				tsipm = ROOT.TMath.Sqrt(float(data[1])**2 - self.sigmatds0[0]**2)
				
				#### Uncertainty contribution from sigma tds0 is negligible. ~ 9.5E-5 ns
				res[fixed_ch]=(float(data[1]), float(data[2]))
		sorted_tuples=sorted(res.items(), key=lambda x:x[1][0])
		sorted_d={k:v for k,v in sorted_tuples}
		self.timeresolutionvalues=sorted_d

	def GetPolyParams(self, runNr, fixed_ch, state='uncorrected', n=5):
		iteration=0 if state=='uncorrected' else 1
		if not os.path.exists(f'{self.path}Polyparams/run{runNr}/polyparams{n}_{fixed_ch}.csv'): return 
		with open(f'{self.path}Polyparams/run{runNr}/polyparams{n}_{fixed_ch}.csv', 'r') as f:
			reader=csv.reader(f)
			alldata=[r for r in reader]
			if len(alldata)==0:return
			data=alldata[iteration]
		
		if n==4:		
			params=[float(i) for i in data[1:11]]
			limits=[float(i) for i in data[11:13]]
			# if len(data)==15: tds0mean=[float(i) for i in data[13:]]
		elif n==5:
			params=[float(i) for i in data[1:13]]
			limits=[float(i) for i in data[13:15]]
			# tds0mean=[float(i) for i in data[15:]]
		return params,limits

	def Gettds0relativetime(self, runNr, fixed_ch, mode='mean', state='uncorrected', n=5):
		
		iteration=0 if state=='uncorrected' else 1
		if not os.path.exists(f'{self.path}Polyparams/run{runNr}/polyparams{n}_{fixed_ch}.csv'): return 
		with open(f'{self.path}Polyparams/run{runNr}/polyparams{n}_{fixed_ch}.csv', 'r') as f:
			reader=csv.reader(f)
			alldata=[r for r in reader]
			if len(alldata)==0:return
			data=alldata[iteration]

		if n==4:
			if len(data)!=15:
				print(f'No alignment constant saved for {fixed_ch}')
				return 
			else: tds0tSiPMmean=[float(i) for i in data[13:]]
		elif n==5:
			if len(data)!=17:
				print(f'No alignment constant saved for {fixed_ch}')
				return 
			else: tds0tSiPMmean=[float(i) for i in data[15:]]
		else: print(f'No tw correction parameters stored for n={n}')
		return tds0tSiPMmean

	def Gettimeresolution(self, runNr, fixed_ch, state, mode='json'):
		
		if mode=='csv':
			fname=f'{self.path}TimeResolution/run{runNr}/timeresolution_{fixed_ch}.csv'
			if not os.path.exists(fname): return 
			with open(fname, 'r') as f:
				reader=csv.reader(f)
				alldata=[row for row in reader]
				if len(alldata)<iteration+1: return 
				data=alldata[iteration]
				# if math.isnan(float(data[0])): return -999.
			timeresolution=float(data[1]), float(data[2])
			return timeresolution

		elif mode=='json':
			filename=f'{self.path}TimeResolution/run{runNr}/timeresolution_{fixed_ch}.json'

			if not os.path.exists(filename): return 

			with open(filename, 'r') as x:
				d=json.load(x)

			if state not in d: return 

			return d[state][0], d[state][1]			

	def FitForMPV(self, runNr, fixed_ch, state):
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
		# xmin,xmax=xproj.GetXaxis().GetXmin(), xproj.GetXaxis().GetXmax()
		mode=xproj.GetBinCenter(xproj.GetMaximumBin())
		res=self.fit_langau(xproj, 'LQ',max(mode-2, 0))
		# 'LQ',0.8*tmp.GetBinCenter(bmin),1.5*tmp.GetBinCenter(bmax)
		MPV=res.Parameter(1)
		MPV_err=res.ParError(1)
		chi2, NDF= res.Chi2(), res.Ndf()
		return (MPV, MPV_err, chi2, NDF)

	def Getchi2_info(self, runNr, fixed_ch, state, n=5):

		fname=f'{self.path}chi2s/run{runNr}/chi2s{n}_{fixed_ch}.csv'
		if not os.path.exists(fname): return 
		with open(fname, 'r') as handle:
			reader=csv.reader(handle)
			alldata=[row for row in reader]
			if len(alldata)==0 or len(alldata)<iteration: return 
			data=alldata[iteration-1]
			try:
				chi2info=str(data[0]), float(data[1]), int(data[2])
			except ValueError:
				print(f'Non-integer NDF for {fixed_ch}')
				return 
		return chi2info

	def Getchi2pNDF(self, runNr, fixed_ch, state, n=5):

		iteration=0 if state=='uncorrected' else 1   
		fname=f'{self.path}chi2s/run{runNr}/chi2s{n}_{fixed_ch}.csv'
		if not os.path.exists(fname): return 
		with open(fname, 'r') as handle:
			reader=csv.reader(handle)
			alldata=[row for row in reader]
			if len(alldata)==0 or len(alldata)<iteration: return 
			data=alldata[iteration-1]
			try:
				chi2pNDF=float(data[1])/int(data[2])
			except ValueError:
				print(f'Non-integer NDF for {fixed_ch}')
				return 
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

	def Makechi2pNDFdict(self, runNr, subsystem, state, n=5):

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

	def GetToFcorrection(self, SiPM, pred, cs, xref):
		c_SiPM=float(cs[0])
		ToFcorrection=abs((pred-xref)/c_SiPM)
		return ToFcorrection

	"""
	Important note: the Analysis.correct_ToF function corrects 
	the SiPM time to the centre of the bar! 
	"""

	def correct_ToF(self, fixed_ch, clock, pred):
		detID=int(fixed_ch.split('_')[0])
		s,p,b=self.parseDetID(detID)
		SiPM=int(fixed_ch.split('_')[-1])		
		xref=self.xrefs[s]
		cs=self.cscintvalues[fixed_ch]

		# fixed_subsystem, fixed_plane, fixed_bar, fixed_SiPM = fixed
		time=clock*self.TDC2ns
		c_SiPM=float(cs[0])
		ToFcorrection=abs((pred-xref)/c_SiPM)
		# Does not work for DS! 
		if SiPM<8:
			if pred >= xref: corrected_t = time - ToFcorrection
			else: corrected_t = time + ToFcorrection
		else: 
			if pred >= xref: corrected_t = time + ToFcorrection
			else: corrected_t = time - ToFcorrection
		return (SiPM, corrected_t)

	def GetCorrectedTimes(self, hit,x=0, mode='unaligned'):
		detID=hit.GetDetectorID()

		alignedtimes=[]
		clocks, qdcs=hit.GetAllTimes(), hit.GetAllSignals()
		for i in clocks:
			SiPM, clock=i 
			fixed_ch=f'{detID}_{SiPM}'
			qdc=self.GetChannelVal(SiPM, qdcs)
			correctedtime=self.MuFilterCorrectedTime(fixed_ch, qdc, clock, x, mode)
			if not correctedtime: continue
			alignedtimes.append((SiPM, correctedtime))
		return alignedtimes

	def MuFilterCorrectedTime(self, fixed_ch, qdc, clock, x=0, mode='unaligned'):
		time=clock*self.TDC2ns

		if not fixed_ch in self.twparameters:
			return

		if not fixed_ch in self.alignmentparameters and mode=='unaligned':
			return
		
		if fixed_ch not in self.cscintvalues and x!=0:
				return
		elif fixed_ch in self.cscintvalues and x!=0: 
			cdata=self.cscintvalues[fixed_ch] 

		twparams=self.twparameters[fixed_ch]
		
		#### Correct ToF if needed.
		if x==0:
			ToFcorrectedtime=time
		else: 
			SiPM=int(fixed_ch.split('_')[-1])
			cscint=self.cscintvalues[fixed_ch]
			xref=self.xrefs[int(fixed_ch[0])]
			ToFcorrectedtime=self.correct_ToF(fixed_ch, clock, x)[1]

		#### TW corrected time then ToF & TW corrected time
		twcorrection=self.correctionfunction(twparams, qdc)
		ToFTWcorrectedtime=ToFcorrectedtime+twcorrection
		
		return ToFTWcorrectedtime

	def GetTimingDiscriminant(self):

		hits=self.task.M.eventTree.Digi_MuFilterHits
		US1hits=[h.GetDetectorID() for h in hits if all([h.GetDetectorID()//10000==2, self.parseDetID(h.GetDetectorID())[1]==0])]

		if len(US1hits)==0: return -999

		elif len(US1hits)==1:
			x=US1hits[0]
			tmp={h.GetDetectorID():h for h in hits}
			us1hit=tmp[x]

		elif len(US1hits)>1:
			docas={}
			for US1detID in US1hits:
				self.task.MuFilter.GetPosition(US1detID, self.task.A, self.task.B)
				docas[self.Getyresidual(US1detID)]=US1detID
			x=docas.pop(min(docas))
			tmp={h.GetDetectorID():h for h in hits}
			us1hit=tmp[x]  
		else:
			print('shite')
			return -999

		averageUS1time=self.GetAverageTime(us1hit)
		if not averageUS1time: return -998

		DS3Haverage=self.GetDSHaverage(hits, mode='timingdiscriminant')

		return DS3Haverage-averageUS1time

	def GetQDCpeak(self, runNr, fixed_ch):
		
		filename=f'{self.path}rootfiles/run{runNr}/timewalk_{fixed_ch}.root'
		if not os.path.exists(filename):return 

		f=ROOT.TFile.Open(filename, 'READ')
		histname=f'dtvqdc_{fixed_ch}_corrected'
		if not hasattr(f, histname): 
			f.Close()
			return 
		hist=f.Get(histname)
		xproj=hist.ProjectionX()
		qdcpeak=xproj.GetBinLowEdge(xproj.GetMaximumBin())
		error=xproj.GetStdDev()/ROOT.TMath.Sqrt(xproj.GetEntries())
		f.Close()
		return qdcpeak, error

	def GetCutDistributions(self, runNr, distmodes=('dy', 'slopes', 'nSiPMs', 'timingdiscriminant'), task='TimeWalk'):
		Allmodes=('dy', 'nSiPMs', 'slopes', 'timingdiscriminant')
		filename=f'{self.path}rootfiles/run{runNr}/SelectionCriteria.root'
		if not os.path.exists(filename): 
			if self.timealignment=='old': filename=f'{self.path}rootfiles/run005097/SelectionCriteria.root'
			else: filename=f'{self.path}rootfiles/run005408/SelectionCriteria.root'

		if isinstance(distmodes, str):
			distmodes=(distmodes,)

		for distmode in distmodes:
			if distmode not in Allmodes:
				print('Use a valid mode.')
				return -999

		f=ROOT.TFile.Open(filename, 'READ')
		dists={}

		for distmode in distmodes:
			if distmode=='dy' or distmode=='nSiPMs':
				for s in (1,2):
					for p in range(self.systemAndPlanes[s]):
						key=str(s*10+p)
						name=f'{distmode}_{key}'
						if not hasattr(f, name):
							print(f'No hist {name}')
							f.Close()
							return -999
						hist=f.Get(name).Clone()
						name=hist.GetName()
						if task=='SelectionCriteria': hist.SetName(f'sc-{name}')
						hist.SetDirectory(ROOT.gROOT)
						dists[name]=hist

			elif distmode=='slopes':
				if not hasattr(f, distmode):
					f.Close()
					print(f'No {distmode} hist')
					return -999
				hist=f.Get(f'{distmode}').Clone()
				hist.SetDirectory(ROOT.gROOT)
				dists[distmode]=hist

			elif distmode=='timingdiscriminant':
				if not hasattr(f, distmode):
					f.Close()
					print(f'No {distmode} hist')
					return -999
				hist=f.Get(f'{distmode}').Clone()
				name=hist.GetName()
				if task=='SelectionCriteria': hist.SetName(f'sc-{name}')
				hist.SetDirectory(ROOT.gROOT)
				dists[distmode]=hist				

		f.Close()
		return dists
	
	def GetTimeAlignmentType(self, runNr):
		if not isinstance(runNr, int): runNr=int(runNr)
		
		if any( [int(runNr) < 5116, int(runNr) > 5174 and int(runNr) < 5193] ): return 'old'
		elif int(runNr) < 5413: return 'new'
		else: return 'new+LHCsynch'

	def isLHCsynch(self, runNr):
		if int(runNr)<5431: return False
		else: return True

	def GetSiPMtime(self, runNr, fixed_ch, state, mode='mean'):
		
		filename=f'{self.path}rootfiles/run{runNr}/timewalk_{fixed_ch}.root'
		if not os.path.exists(filename): return 

		f=ROOT.TFile.Open(filename, 'READ')
		histname=f'tSiPM_{fixed_ch}_{state}'
		if not hasattr(f, histname): 
			f.Close()
			return 
		hist=f.Get(histname)
		if mode=='mean':
			tSiPM=hist.GetMean()
			uncertainty=hist.GetStdDev()
		elif mode=='median':
			pass
		f.Close()
		return tSiPM, uncertainty

	def GetAlignmentParameters(self, runNr, fixed_ch):
		fname=f'{self.path}Alignmentparams/run{runNr}/alignmentparams{self.CorrectionType}_{fixed_ch}.csv'
		if not os.path.exists(fname): 
			return
		with open(fname, 'r') as handle:
			reader=csv.reader(handle)
			alldata=[row for row in reader]
			if len(alldata) == 0: return 
			try:
				data=alldata[0]
			except IndexError:
				print(f'{fixed_ch} IndexError')
		return (float(data[0]), float(data[1]))		

	def Gettds0mean(self, runNr, fixed_ch, mode='mean', state='corrected'):
		filename=f'{self.path}rootfiles/run{runNr}/timewalk_{fixed_ch}.root'
		if not os.path.exists(filename): 
			print(f'No timewalk file for {fixed_ch}')
			return

		f=ROOT.TFile.Open(filename, 'READ')
		# if state=='uncorrected': histname=f'dtvqdc_{fixed_ch}_uncorrected'
		# elif state=='corrected': histname=f'dtvqdc_{fixed_ch}_corrected'
		if state=='aligned': histname=f'AlignedSiPMtime_{fixed_ch}'
		else: print(f'No conditions met for finding histname:\n')

		if not hasattr(f, histname): 
			f.Close()
			print(f'No hist {histname} for {fixed_ch}')
			return 
		hist=f.Get(histname)

		if mode=='mean':
			mean=hist.GetMean()
			uncertainty = mean/ROOT.TMath.Sqrt(hist.GetEntries())
			f.Close()
			return mean, uncertainty

		elif mode=='truncated':
			modaltime=hist.GetBinCenter(hist.GetMaximumBin())
			low, high=hist.FindBin(modaltime-hist.GetStdDev()), hist.FindBin(modaltime+hist.GetStdDev()) # mode(hist) + 2*stddev, mode(hist) - 2*std dev
			tmp=[[hist.GetBinLowEdge(i),hist.GetBinContent(i)] for i in range(low, high+1)] # Data within 2 standard deviations of the 
			mean=sum([x[0]*x[1] for x in tmp]) / sum([x[1] for x in tmp])
			# uncertainty = mean / ROOT.TMath.Sqrt(sum( [x[1] for x in tmp] ))
			uncertainty = ROOT.TMath.Sqrt(1/len(tmp) * sum([(x[0]-mean)**2 for x in tmp]))
			f.Close()
			return mean, uncertainty      

		elif mode=='mode':
			correction=hist.GetBinCenter(hist.GetMaximumBin())
			uncertainty=hist.GetStdDev()/hist.GetEntries()
			return correction, uncertainty
		f.Close()

	def GettSiPMcorrectedmean(self, runNr, fixed_ch, mode='mean'):
		filename=f'{self.path}rootfiles/run{runNr}/timewalk_{fixed_ch}.root'
		if not os.path.exists(filename): return
		f=ROOT.TFile.Open(filename, 'READ')
		# if self.timealignment=='new': histname=f'tSiPMToFcorrected_{fixed_ch}_corrected_tDS0-tSiPMcorrected'
		histname=f'tSiPMToFcorrected_{fixed_ch}_corrected_tDS0-tSiPMcorrected'
		# elif self.timealignment=='old': histname=f'ScifiAlignedToFcorrectedtSiPM_{fixed_ch}_corrected'
		if not hasattr(f, histname): 
			f.Close()
			return 
		hist=f.Get(histname)
		if mode=='mean':
			correction=hist.GetBinCenter(hist.GetMaximumBin())
			error=hist.GetStdDev()/2
			# for 
			# correction=hist.GetMean()

			# uncertainty=hist.GetStdDev()

		elif mode=='median':
			correction, error =0, 0
		else: 
			correction, error =0, 0
		f.Close()
		return correction, error

	def GetEntriesInHist(self, runNr, fixed_ch, mode, state):

		filename=f'{self.path}rootfiles/run{runNr}/timewalk_{fixed_ch}.root'
		if not os.path.exists(filename):return 

		f=ROOT.TFile.Open(filename, 'READ')
		histname=f'{mode}_{fixed_ch}_{state}'
		if not hasattr(f, histname): 
			f.Close()
			return 
		hist=f.Get(histname)
		entries=hist.GetEntries()
		f.Close()
		return entries

	# def GetPolyParamRanges(self, runNr, subsystem, state, n=5):
	# 	params={i:[0.,0.] for i in ('A', 'B', 'C')}
	# 	path=f'{self.path}Polyparams/run{runNr}/'
	# 	files=[i for i in os.listdir(path) if int(i.split('_')[1][0])==subsystem]
	# 	counter=0
	# 	for f in files:
	# 		if f.split('_')[0].find(str(n)) == -1:continue
	# 		#if idx>1: break
	# 		fixed_ch=f"{f.split('_')[1]}_{f.split('.')[0].split('_')[-1]}"
	# 		tmp = self.GetPolyParams(runNr, fixed_ch, n, iteration)
	# 		if isinstance(tmp, int): continue
	# 		else: fps=tmp[0]
	# 		if fps==-999.: 
	# 			print(f'{fixed_ch} has no poly params')
	# 			continue

	# 		if n==1:
	# 			vals={p:fps[i*2] for i,p in enumerate(('A', 'B', 'C'))}
	# 			if counter==0:
	# 				for i in ('A', 'B', 'C'): params[i]=[vals[i], 1.1*vals[i]]
	# 				# for idx,x in enumerate(('A', 'B', 'C')): params[x]=[fps[idx*2], fps[idx*2+1]]
	# 		else: 
	# 			for i in ('A', 'B', 'C'):
	# 				if vals[i] < params[i][0]: params[i][0]=vals[i]
	# 				if vals[i] > params[i][1]: params[i][1]=1.1*vals[i]
	# 			if n==5:
	# 			vals={p:fps[i*2] for i,p in enumerate(('A', 'B', 'C'))}
	# 			if counter==0:
	# 				for i in ('A', 'B', 'C'): params[i]=[vals[i], 1.1*vals[i]]
	# 				# for idx,x in enumerate(('A', 'B', 'C')): params[x]=[fps[idx*2], fps[idx*2+1]]
	# 			else: 
	# 				for i in ('A', 'B', 'C'):
	# 					if vals[i] < params[i][0]: params[i][0]=vals[i]
	# 					if vals[i] > params[i][1]: params[i][1]=1.1*vals[i]         
	# 		counter+=1

	# 	return params

	def MakeTWCorrectionDict(self, runNr, withErrors=False):
  
		d={}
		for s in (1,2):
			for p in range(self.systemAndPlanes[s]):
				for b in range(self.systemAndBars[s]):
					for SiPM in self.systemAndSiPMs[s]:
						fixed_ch=self.MakeFixedCh((s,p,b,SiPM))
						tmp=self.GetPolyParams(runNr, fixed_ch, state='uncorrected', n=self.CorrectionType)
						if not tmp: continue
						else: paramsAndErrors=tmp[0]
						if not withErrors: params=self.correctionparams(paramsAndErrors)
						else: params=paramsAndErrors
						d[fixed_ch]=params
		if len(d)==0: self.twparameters=None
		self.twparameters=d

	### Make dictionary of the alignment parameter determined as the truncated y-mean of tw-corr (tds0 - tSiPM)
	def MakeAlignmentParameterDict(self, runNr):
		d={}
		for s in (1,2):
			for p in range(self.systemAndPlanes[s]):
				for b in range(self.systemAndBars[s]):
					for SiPM in self.systemAndSiPMs[s]:
						fixed_ch=self.MakeFixedCh((s,p,b,SiPM))
						# correction=self.Gettds0mean(runNr, fixed_ch)
						correction=self.GetAlignmentParameters(runNr, fixed_ch)
						if not correction: continue
						d[fixed_ch]=correction
		if len(d)==0: self.alignmentparameters=None
		self.alignmentparameters=d

	def MakeTimingCovarianceDict(self, runNr):
		filename=f'{self.path}TimingCovariance/run{self.runNr}/timingcovariance.json'
		with open(filename) as jsonfile:
			self.timingcovariance=json.load(jsonfile)

	def MakeTimingCorrelationDict(self, runNr):
		filename=f'{self.path}TimingCovariance/run{self.runNr}/timingcorrelation.json'
		with open(filename) as jsonfile:
			self.timingcorrelation=json.load(jsonfile)

	def MakeQDCMIPJson(self, runNr):
		### Currently 
		self.langaufun()
		d={}
		for s in (1,2):
			for p in range(self.systemAndPlanes[s]):
				for b in range(self.systemAndBars[s]):
					print(f'{s}, {p}, {b}')
					for SiPM in self.systemAndSiPMs[s]:
						fixed_ch=self.MakeFixedCh((s,p,b,SiPM))
						MPV, MPV_err, chi2, NDF = self.FitForMPV(runNr, fixed_ch, 'corrected')

						d[fixed_ch] = MPV, MPV_err

		mpvpath=f'{self.afswork}MPVs/run{self.runNr}/'
		mpvfilename=covariancepath+f'MPVs.json'
		if not os.path.exists(mpvfilename): 
			os.makedirs(mpvpath, exist_ok=True)
		with open(mpvfilename, 'w') as outfile:
			json.dump(d, outfile, indent=4)
		print(f'mpv dictionary written to {mpvfilename}')

	def GetAverageSiPMTimingCovariance(self, runNr):
		r={}
		
		if not hasattr(self, 'timingcovariance'): self.MakeTimingCovarianceDict(runNr)

		for key in self.timingcovariance:
			if key=='timingxt_tds0_Veto' or key =='timingxt_tds0_US': continue
			x, detID, tmp=key.split('_')
			
			SiPMs=tmp[len('SiPMs'):].split('-')

			#### Special case for last SiPM on bar end due to how combinations algorithm works
			if SiPMs[1] in ('7', '15'):
				fixed_ch=f'{detID}_{SiPMs[1]}'
				if not fixed_ch in r: r[fixed_ch]=[]
				r[fixed_ch].append(self.timingcovariance[key])

			fixed_ch = f'{detID}_{SiPMs[0]}'
			if not fixed_ch in r: r[fixed_ch]=[]
			r[fixed_ch].append(self.timingcovariance[key])

		for fixed_ch in r: 
			avg=1/len(r[fixed_ch]) * sum(r[fixed_ch]) 
			r[fixed_ch]=avg 
		
		self.SiPMaveragetimingcovariance=r
	
	def GetCovariance(self, runNr, detID, SiPMs):
		if not hasattr(self, 'timingcovariance'):
			self.MakeTimingCovarianceDict(runNr)

		key=f'timingxt_{detID}_SiPMs{SiPMs[0]}-{SiPMs[1]}'
		if not key in self.timingcovariance: return
		else: return self.timingcovariance[key]

	def GetCanvas(self, runNr, fixed_ch, mode, iteration):
		filename=f'{self.path}rootfiles/run{runNr}/timewalk_{fixed_ch}.root'
		f=ROOT.TFile.Open(filename, 'read')
		name=f'{mode}_{fixed_ch}_{iteration}'
		if not hasattr(infile, name):
			print(f'No canvas with name: {name} available in file')
			f.Close()
			return 
		og_canv=f.Get(name)
		canv.og_canv.Clone()
		canv.SetDirectory(ROOT.gROOT)
		f.Close()
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
