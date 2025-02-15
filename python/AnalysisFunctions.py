#!/usr/bin/env python
import ROOT,os,csv,json
from datetime import datetime
import math as m 
import numpy as np
from itertools import combinations

# class mufi_analysis: # I could make several classes for generic functions like parseDetID and more specific ones for 

class Analysis(object):

	def __init__(self, options):

		self.options=options

		# Adding flag for LaserMeasurements/ work to use some functions defined in here
		if not hasattr(options, "LaserMeasurements"):
			if options.runNumber==-1: self.runNr='005408'
			else: self.runNr = str(options.runNumber).zfill(6)
			self.TWCorrectionRun = str(5408).zfill(6)

			self.timealignment=self.GetTimeAlignmentType(self.runNr)
			self.state=options.state
			if hasattr(options, 'datafiletype'): self.fileext=options.datafiletype
			else: self.fileext='json'		
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

			self.referencesystem=options.referencesystem
			self.refsysname='DS' if self.referencesystem==3 else 'SF'
			
			if options.numuStudy: self.Get_numuevents()	
			elif options.nueStudy: self.Get_nueevents()		
			
			self.simulation = options.simulation
			if self.simulation: 
				self.simEngine = self.GetSimEngine

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
		freq = 160.316E6
		self.TDC2ns = 1E9/freq
			
	def GetSimEngine(self):
		simEngine = self.options.geoFile.split('.')[1].split('-')[0]
		return simEngine 
	
	def print_timestamp(self, message=""):
		print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} - {message}")

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
		elif 5485 <= runNr < 7656 : year='2023'
		else: year='2024'

		return year

	def BuildBarLengths(self, MuFilter):
		Vetobarlength = MuFilter.GetConfParF('MuFilter/VetoBarX')
		USbarlength = MuFilter.GetConfParF('MuFilter/UpstreamBarX')
		DSbarlength_hor = MuFilter.GetConfParF('MuFilter/DownstreamBarX')
		DSbarlength_vert = MuFilter.GetConfParF('MuFilter/UpstreamBarY_ver')
		
		self.barlengths={1:Vetobarlength, 2:USbarlength, 3:DSbarlength_hor, 4:DSbarlength_vert}
		
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

	def GetScifiAverageTime(self, scifi, scifihits, mode='tsf0'):

		stations = {i.GetDetectorID():i for i in scifihits}
		times=[]	

		if mode=='tsf0':
			for hit in scifihits:
				detID = hit.GetDetectorID()
				t = scifi.GetCorrectedTime(detID, hit.GetTime(), 0)
				times.append(t)
			if len(times)==0: return False 
			else: return sum(times) / len(times)
		
		elif mode=='deltastations':
			nstations = scifi.GetConfParI('Scifi/nscifi')
			times={st:[] for st in range(nstations+1)}
			res={f'scifi-delta{i+1}{i}':0 for i in range(nstations)}
			for hit in scifihits:
				detID = hit.GetDetectorID()
				station = hit.GetStation()
				t = scifi.GetCorrectedTime(detID, hit.GetTime(), 0)
				times[station].append(t)
			for station in range(nstations-1):
				if len(times[station+1]) == 0 or len(times[station]) == 0: continue
				res[f'scifi-delta{station+1}{station}'] = sum(times[station+1])/len(times[station+1]) - sum(times[station])/len(times[station])
			return res
			
	def GetScifiTrackAverageTime(self, scifi, scifihits):

		stations = {i.GetDetectorID():i for i in scifihits}
		times=[]

		track = self.task.track
		nM = track.getNumPointsWithMeasurement()

		for n in range(nM):
			M=track.getPointWithMeasurement(n)
			W=M.getRawMeasurement()
			detID=W.getDetId()
	
			# hkey=W.getHitId()
			if detID not in stations:
				print(stations)
				print(f'Event {self.task.M.EventNumber}')
			trackHit = stations[detID]
			time = scifi.GetCorrectedTime(detID, trackHit.GetTime(), 0)
			times.append(time)
		
		if len(times)==0: return 
		averagetime = sum(times) / len(times)
		return averagetime

	def GetPlaneData(self, hits, mufilter):
		planewise_data = {i:{} for i in range(5)}

		# Group hits by plane
		for hit in hits:

			detID = hit.GetDetectorID()
			s,p,b = self.parseDetID(detID)
			if not s==2: continue
			if not hit.isValid(): continue

			# muAna knows the reference time needed to apply for !simulation
			if not self.simulation: alignedtimes=self.GetCorrectedTimes(hit, mode='aligned')
			else: alignedtimes=hit.GetAllTimes()

			atimes_left = [i[1] for i in alignedtimes if i[0]<8]
			atimes_right = [i[1] for i in alignedtimes if i[0]>=8] 

			# Check nSiPMs left and right. Skip hit if less than 4 SiPMs on either left or right fire
			if len(atimes_left)<4 or len(atimes_right)<4:continue

			atimes_left_mean = sum(atimes_left)/len(atimes_left)
			atimes_right_mean = sum(atimes_right)/len(atimes_right)

			# Skip hit if abs( atimes_left(right) ) > L/2 / cscint_L(R)
			averagecscint_left, averagecscint_right = self.GetBarAveragecscint(mufilter, detID)
			if abs(atimes_left_mean) > self.barlengths[2]/2 / averagecscint_left[0] or abs(atimes_right_mean) > self.barlengths[2]/2 / averagecscint_right[0]: continue

			planewise_data[p][detID] = {}
			planewise_data[p][detID]['bar-QDC'] = self.GetTotalQDC(hit.GetAllSignals())

			planewise_data[p][detID]['atimes-left'] = atimes_left_mean
			planewise_data[p][detID]['atimes-right'] = atimes_right_mean

			planewise_data[p][detID]['cscint-left'] = averagecscint_left
			planewise_data[p][detID]['cscint-right'] = averagecscint_right	

		return planewise_data	

	def GetBarycentres(self, hits, **kwargs):

		"""
		Adding some kwargs for using when the analysis instance
		isn't connected to a FairTask
		"""

		mufilter=kwargs.get("MuFilter")
		if "hasTrack" in kwargs: hasTrack=kwargs.get("hasTrack")
		else: hasTrack=False

		if not hasattr(self, "barlengths"): self.BuildBarLengths(mufilter)

		barycentres={i:{} for i in range(5)}

		# Plane data rejects invalid hits! 
		planewise_data = self.GetPlaneData(hits, mufilter)

		for plane in planewise_data:
			pdata = planewise_data[plane]
			if len(pdata)==0: continue

			# Determine quantities for the bars now
			planeQDC = sum([ pdata[detID]['bar-QDC'] for detID in pdata ])

			weighted_ys = []
			x_barycentres = {}

			for detID in pdata:
				# Get weighted y-position for each hit that passes selection in this plane
				s,p,b = self.parseDetID(detID)
				
				if hasTrack:
					if self.GetExtrapolatedBarDetID(p) == detID: trackInBar = True
					else: trackInBar=False
				else: trackInBar=False

				barQDC = pdata[detID]['bar-QDC']
				mufilter.GetPosition(detID, self.A, self.B)
				y_pos = barQDC/planeQDC * 0.5 * (self.A.y() + self.B.y())
				x_midpoint = 0.5 * (self.A.x() + self.B.x())
				weighted_ys.append(y_pos) # Get weighted y-pos

				# Get x-barycentre
				cscint_left, cscint_right = pdata[detID]['cscint-left'], pdata[detID]['cscint-right']
				atimes_left, atimes_right = pdata[detID]['atimes-left'], pdata[detID]['atimes-right']

				dxLphys, dxRphys = (x_midpoint+atimes_left*cscint_left[0]), (-atimes_right*cscint_right[0]+x_midpoint) # Adding the x_midpoint translates into the physics FoR
				dxLphys_err, dxRphys_err = self.Getxuncertainty(detID, pdata[detID],'left'),self.Getxuncertainty(detID, pdata[detID],'right')

				x_barycentre = 0.5*(dxLphys + dxRphys)
				lambda_x = (dxLphys - dxRphys)

				if trackInBar:

					xEx, yEy, zEx = self.GetExtrapolatedPosition(p)
					d_left_track = dxLphys - xEx
					d_right_track = dxRphys - xEx					

					shower_side_array = np.where(abs(d_left_track) < abs(d_right_track), 'right', 'left')
					
					if len(np.unique(shower_side_array)) == 1:
						shower_side = str(np.unique(shower_side_array)[0])  # Extract the single unique value
					else:
						# Handle the case where shower_side_array has multiple values
						raise ValueError("shower_side_array contains multiple values, unable to determine a single key")
				else: shower_side=False					
				
				# Store barycentre determined by each bar
				x_barycentres[detID] = {
					'xL':(dxLphys,dxLphys_err),
					'xR':(dxRphys,dxRphys_err),
					'xB':x_barycentre,
					'lambda_x':lambda_x,
					'barQDC':barQDC,
					"trackInBar":shower_side
					}

			y_barycentre=sum(weighted_ys) 
				
			mufilter.GetPosition(max(pdata.keys()), self.A, self.B)
			max_y = 0.5*(self.A.y() + self.B.y())
			mufilter.GetPosition(min(pdata.keys()), self.A, self.B)
			min_y = 0.5*(self.A.y() + self.B.y())
			lambda_y = max_y - min_y
			
			y_barycentres = {
				'yB':y_barycentre,
				'lambda_y':lambda_y
				}

			barycentres[plane] = {'x-barycentres':x_barycentres, "y-barycentre":y_barycentres}				

		return barycentres
	
	def GetOverallXBarycentre(self, barycentres, mode):
		xs = {p:{} for p in barycentres.keys()} 
		
		if mode=='relQDC':
			for p,pdata in barycentres.items():
				if len(pdata)==0: continue
				xb_data = pdata['x-barycentres']
				
				planeQDC = sum([xb_data[detID]['barQDC'] for detID in xb_data])
				relQDCs={detID:xb_data[detID]['barQDC']/planeQDC for detID in xb_data.keys()}
				
				xL = sum([relQDCs[detID]*xb_data[detID]['xL'][0] for detID in xb_data.keys()])
				sigma_xL = np.sqrt(sum([relQDCs[detID]**2*xb_data[detID]['xL'][1]**2 for detID in xb_data.keys()]))

				xR = sum([relQDCs[detID]*xb_data[detID]['xR'][0] for detID in xb_data.keys()])
				sigma_xR = np.sqrt(sum([relQDCs[detID]**2*xb_data[detID]['xR'][1]**2 for detID in xb_data.keys()]))

				xB = sum([relQDCs[detID]*xb_data[detID]['xB'] for detID in xb_data.keys()])
				
				lambda_x = sum([relQDCs[detID]*xb_data[detID]['lambda_x'] for detID in xb_data.keys()])

				xs[p]['dxL']=(xL, sigma_xL)
				xs[p]['dxR']=(xR, sigma_xR)
				xs[p]['dxB']=xB
				xs[p]['lambda_x']=lambda_x
			return xs
		
		elif mode=='maxQDC':
			for p,pdata in barycentres.items():
				if len(pdata)==0: continue
				xb_data=pdata['x-barycentres']
				
				max_key = max(xb_data, key=lambda k: xb_data[k]['barQDC'])

				xL=xb_data[max_key]['xL']
				xR=xb_data[max_key]['xR']
				xB=xb_data[max_key]['xB']
				lambda_x=xb_data[max_key]['lambda_x']

				xs[p]['dxL']=xL
				xs[p]['dxR']=xR
				xs[p]['dxB']=xB
				xs[p]['lambda_x']=lambda_x
			return xs
	
	def Getxuncertainty(self, detID, bardata, side):

		barside_sigma=self.CalculateBarsideTimeresolution(self.runNr, detID, side)

		ct=bardata[f'cscint-{side}'][0]*bardata[f'atimes-{side}']
		c_rel_error = bardata[f'cscint-{side}'][1]/bardata[f'cscint-{side}'][0]
		
		# print(bardata[showerside])
		t_rel_error = barside_sigma/bardata[f'atimes-{side}']

		dx_sq = ct**2 * ( c_rel_error**2 + t_rel_error**2 )
		dx = np.sqrt(dx_sq)
		return dx
		
	def GetDSHaverage(self, hits, mode='tds0'):

		"""
		For tds0, I am only using DS2 and 3 due to the delta tds21 plot being broad and asymmetric... not implemented, should I?
		"""
		stations = {i.GetDetectorID():i for i in hits}

		if mode=='tds0':
			total={i:0 for i in range(4)}
			counter={i:0 for i in range(4)}
			theTrack = self.task.M.Reco_MuonTracks[0]
			nM = theTrack.getNumPointsWithMeasurement()
			# print(f'{nM} measurements')
			for n in range(nM):
				M=theTrack.getPointWithMeasurement(n)
				W=M.getRawMeasurement()
				detID=W.getDetId()
				s,p,b = self.parseDetID(detID)
				hkey=W.getHitId()
				trackHit = stations[detID]
				tdcs = trackHit.GetAllTimes()
				# print(f'DetID: {detID}, len(tdcs)={len(tdcs)}')
				if not all ( [p in (0,1,2), b<60, len(tdcs)==2] ): continue

				if len(tdcs) != 2: continue # pass

				for item in tdcs: 
					SiPM, clock = item
					dscorrectedtime=self.task.MuFilter.GetCorrectedTime(detID, SiPM, clock*self.TDC2ns, 0)
					total[p]+=dscorrectedtime
					counter[p]+=1

			# Only using DS3 and DS2
			sum_total, counter_total = total[2] + total[1], counter[2] + counter[1]
			# print(f'total, counter: {sum_total}, {counter_total}')
			if counter_total==0: return -999

			dsh_average, nfired = sum_total/counter_total, counter_total

			return dsh_average, nfired	

		elif mode=='timingdiscriminant': # Take k=2 for the 3rd horizontal plane

			total=[]

			theTrack = self.task.M.Reco_MuonTracks[0]
			nM = theTrack.getNumPointsWithMeasurement()
			# print(f'{nM} measurements')
			for n in range(nM):
				M=theTrack.getPointWithMeasurement(n)
				W=M.getRawMeasurement()
				detID=W.getDetId()
				s,p,b = self.parseDetID(detID)
				hkey=W.getHitId()
				trackHit = stations[detID]
				tdcs = trackHit.GetAllTimes()
				# print(f'DetID: {detID}, len(tdcs)={len(tdcs)}')
				if not all ( [p in (0,1,2), b<60, len(tdcs)==2] ): continue

				if len(tdcs) != 2: continue # pass

				for item in tdcs: 
					SiPM, clock = item
					dscorrectedtime=self.task.MuFilter.GetCorrectedTime(detID, SiPM, clock*self.TDC2ns, 0)
					total.append(dscorrectedtime)

			if len(total)==0:
				return -420					
			return sum(total)/len(total)

		elif mode=='deltastations':
			res={f'delta32':0, 'delta21':0}
			ts={k:[] for k in range(3)}

			theTrack = self.task.M.Reco_MuonTracks[0]
			nM = theTrack.getNumPointsWithMeasurement()
			# print(f'{nM} measurements')
			for n in range(nM):
				M=theTrack.getPointWithMeasurement(n)
				W=M.getRawMeasurement()
				detID=W.getDetId()
				s,p,b = self.parseDetID(detID)
				hkey=W.getHitId()
				trackHit = stations[detID]
				tdcs = trackHit.GetAllTimes()
				# print(f'DetID: {detID}, len(tdcs)={len(tdcs)}')
				if not all ( [p in (0,1,2), b<60, len(tdcs)==2] ): continue

				if len(tdcs) != 2: continue # pass

				for item in tdcs: 
					SiPM, clock = item
					dscorrectedtime=self.task.MuFilter.GetCorrectedTime(detID, SiPM, clock*self.TDC2ns, 0)
					ts[p].append(dscorrectedtime)

			if any([len(ts[i])==0 for i in ts]): return False
			res['delta32']=sum(ts[2])/len(ts[2]) - sum(ts[1])/len(ts[1])
			res['delta21']=sum(ts[1])/len(ts[1]) - sum(ts[0])/len(ts[0])
			return res

		elif mode=='testing-tds0': 
			total={i:0 for i in range(4)}
			counter={i:0 for i in range(4)}
			theTrack = self.task.M.Reco_MuonTracks[0]
			for nM in range(theTrack.getNumPointsWithMeasurement()):
				M=theTrack.getPointWithMeasurement(nM)
				W=M.getRawMeasurement()
				detID=W.getDetId()
				s,p,b = self.parseDetID(detID)
				hkey=W.getHitId()
				# trackHit=hits[hkey]
				trackHit = stations[detID]
				tdcs = trackHit.GetAllTimes()
				if not all ( [p in (0,1,2), b<60, len(tdcs)!=2] ): continue

				if len(tdcs) != 2: continue # pass

				for item in tdcs: 
					SiPM, clock = item
					dscorrectedtime=self.task.MuFilter.GetCorrectedTime(detID, SiPM, clock*self.TDC2ns, 0)
					total[p]+=dscorrectedtime
					counter[p]+=1

			# Only using DS3 and DS2
			sum_total, counter_total = total[2] + total[1], counter[2] + counter[1]
			if counter_total==0: return -999

			dsh_average, nfired = sum_total/counter_total, counter_total

			return dsh_average, nfired				

		elif mode=='testing-deltastations': 
			res={f'delta32':0, 'delta21':0}
			ts={k:[] for k in range(3)}
			theTrack = self.task.M.Reco_MuonTracks[0]
			for nM in range(theTrack.getNumPointsWithMeasurement()):
				M=theTrack.getPointWithMeasurement(nM)
				W=M.getRawMeasurement()
				detID=W.getDetId()
				s,p,b = self.parseDetID(detID)
				if not all ( [p in (0,1,2), b<60] ): continue

				hkey=W.getHitId()
				trackHit = stations[detID]
				times=trackHit.GetAllTimes()
				if not len(times)==2: continue

				for idx in times:
					SiPM, clock=idx
					dscorrectedtime=self.task.MuFilter.GetCorrectedTime(detID, SiPM, clock*self.TDC2ns, 0)
					ts[p].append(dscorrectedtime)

			if any([len(ts[i])==0 for i in ts]): 
				# print(f'ts: {ts}')
				return False
			res['delta32']=sum(ts[2])/len(ts[2]) - sum(ts[1])/len(ts[1])
			res['delta21']=sum(ts[1])/len(ts[1]) - sum(ts[0])/len(ts[0])
			return res

	def fit_langau(self, hist,o,bmin,bmax):
		params = {0:'Width(scale)',1:'mostProbable',2:'norm',3:'sigma'}
		F = ROOT.TF1('langau',self.langaufun,0,200,4)
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

		return (par[2] * step * summe * invsq2pi / par[3])

	def DSAcceptanceRange(self, MuFilter):
		barNumbers={
			'horizontal':
				{'top':'059', 'bottom':'000', 'left':'029', 'right':'029'},
			'vertical':
				{'top':'089', 'bottom':'089', 'left':'119', 'right':'060'}
			}

		vals = {'horizontal':{}, 'vertical':{}}
		res={}

		# x for horizontal planes, y for vertical planes
		for i in barNumbers.keys():
			for plane in range(4):

				if i=='horizontal' and plane==3:continue # only 4th vertical plane
				
				for position in barNumbers[i]:
					bar = barNumbers[i][position]
					detID=f'3{plane}{bar}'
				
					MuFilter.GetPosition(int(detID), self.A, self.B)
					
					if i=='horizontal' and position in ['top', 'bottom']: 
						vals[i][position] = 1/2 * (self.A.y()+self.B.y())
					elif i=='horizontal' and position == 'left': vals[i][position] = self.A.x()
					elif i=='horizontal' and position == 'right': vals[i][position] = self.B.x()
					
					elif i=='vertical' and position in ['left', 'right']: 
						vals[i][position] = 1/2 * (self.A.x()+self.B.x())
					elif i=='vertical' and position == 'top': vals[i][position] = self.A.y()
					elif i=='vertical' and position == 'bottom': vals[i][position] = self.A.y()

		return vals

	def GetAverageTime(self, mufiHit, side='both', correctTW=True):
		value=[0, 0]
		count=[0, 0]
		nSiPMs=mufiHit.GetnSiPMs()
		if correctTW: times = self.GetCorrectedTimes(mufiHit) # in nanoseconds
		else: times = mufiHit.GetAllTimes() # in clockcycles
		
		detID=mufiHit.GetDetectorID()
		s, p, b = self.parseDetID(detID)
		for element in times:
			SiPM, clock = element

			# Check if this channel has determined alignment parameters
			fixed_ch = f'{detID}_{SiPM}'
			if correctTW and not fixed_ch in self.alignmentparameters: continue 
			elif correctTW and fixed_ch in self.alignmentparameters: d = self.alignmentparameters[fixed_ch]
			if s==2 and self.IsSmallSiPMchannel(SiPM):continue
			
			if SiPM<nSiPMs:
				if not correctTW: value[0]+=clock*self.TDC2ns # in nanoseconds. 
				else: value[0]+= clock - d[0] # in nanoseconds ## I can't subtract d[0] here because then I have a sign error 
				# else: value[0]+= clock # in nanoseconds    
				count[0]+=1
			else:
				if not correctTW: value[0]+=clock*self.TDC2ns
				else: value[0]+= clock - d[0] # in nanoseconds ## I can't subtract d[0] here because then I have a sign error 
				# else: value[0]+= clock
				count[1]+=1
    
		if s == 2:
			if count[0] != 0 and count[1] != 0:
				if side == 'both':
					# average = 0.5*(value[0]/count[0]+value[1]/count[1])
					average = sum(value) / sum(count) 
					return average
				if side == 'L':
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

	def GetTotalQDC(self, signals):
		return sum( [i[1] for i in signals if not self.IsSmallSiPMchannel(i[0])] )

	def GetChannelVal(self, SiPM, chs):
		for entry in chs:
			fSiPM, val = entry
			if fSiPM == SiPM:
				return val
		return
		
	def Get_numuevents(self):
		numusignalevent_filepath = '/afs/cern.ch/work/a/aconsnd/numusignalevents.csv'
		self.nu_mu_events = {}
		
		with open(numusignalevent_filepath, 'r') as f:
			reader=csv.reader(f)
			# next(reader) # skip first row with the headers

			for idx,x in enumerate(reader):
				if idx==0: continue

				self.nu_mu_events[int(x[0])] = [int(x[1]), int(x[2])] + [float(i) for i in x[3:]]
	
	def GetNeutrinoIntType(self, event):

		if not hasattr(event, "MCTrack"):
			print(f'No MCTrack branch. Is this real data?')
			return 


		if event.MCTrack[0].GetPdgCode() == event.MCTrack[1].GetPdgCode():
			i_flav = 0 #NC
		elif abs(event.MCTrack[1].GetPdgCode()) == 11:
			i_flav = 1 #nueCC
		elif abs(event.MCTrack[1].GetPdgCode()) == 13:
			i_flav = 2 #numuCC
		elif abs(event.MCTrack[1].GetPdgCode()) == 15:
			is1Mu = False
			for j_track in range(2, len(event.MCTrack)):
				if event.MCTrack[j_track].GetMotherId() == 1 and abs(event.MCTrack[j_track].GetPdgCode()) == 13:
					is1Mu = True
					break
			if is1Mu:
				i_flav = 4 #nutauCC1mu
			else:
				i_flav = 3 #nutauCC0mu    

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
			geoobject.GetPosition(10000, self.A, self.B)
			xmin, xmax = self.B.x(), self.A.x()
			ymin = self.A.y() - geoobject.GetConfParF('MuFilter/VetoBarY')/2.
			geoobject.GetPosition(10006, self.A, self.B)
			ymax = self.A.y() + geoobject.GetConfParF('MuFilter/VetoBarY')/2.
		elif subsystem==2:
			geoobject.GetPosition(20000, self.A, self.B)
			xmin, xmax = self.B.x(), self.A.x()
			ymin = self.A.y() - geoobject.GetConfParF('MuFilter/UpstreamBarY')/2.
			geoobject.GetPosition(20009, self.A, self.B)
			ymax = self.A.y() + geoobject.GetConfParF('MuFilter/UpstreamBarY')/2.
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

	def AllLiveSiPMs(self, hit):
		detID = hit.GetDetectorID()
		if detID != 21000: left_max=right_max=6
		else: left_max, right_max = 6, 5

		nleft, nright = self.GetnFiredSiPMs(hit)
		if nleft==left_max and nright==right_max: return True 
		else: return False

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

	def GetExtrapolatedBarDetID(self, plane):
		xEx, yEx, zEx = self.GetExtrapolatedPosition(plane)
		detIDEx = self.task.nav.FindNode(xEx, yEx, zEx).GetName()
		if not detIDEx.split('_')[0] == 'volMuUpstreamBar': 
			return
		else: return int(detIDEx.split('_')[1])

	def GetExtrapolatedPosition(self, plane):
		zEx = self.zPos['MuFilter'][20+plane]
		lam = (zEx-self.task.pos.z())/self.task.mom.z()
		yEx = self.task.pos.y() + lam*self.task.mom.y()
		xEx = self.task.pos.x() + lam*self.task.mom.x()
		
		return (xEx, yEx, zEx)

	def Getcscint(self, runNr, fixed_ch, state):

		filename=f'{self.path}cscintvalues/run{runNr}/cscint_{fixed_ch}.json'

		if not os.path.exists(filename): return 

		with open(filename, 'r') as x:
			d=json.load(x)

		if state not in d: return 

		if len(d[state])==6: return d[state][0], d[state][1]
		elif len(d[state])==7: return d[state][0], d[state][1], d[state][-1]

	def Getcscint_offset(self, runNr, fixed_ch, state):

		filename=f'{self.path}cscintvalues/run{runNr}/cscint_{fixed_ch}.json'

		if not os.path.exists(filename): return 

		with open(filename, 'r') as x:
			d=json.load(x)

		if state not in d: return 

		return d[state][2], d[state][3]

	def GetBarAveragecscint(self, MuFilter, detID):

		subsystem, plane, bar = self.parseDetID(detID)
		cscintvalues={'left':[], 'right':[]}
		for SiPM in self.systemAndSiPMs[subsystem]:
			fixed_ch=self.MakeFixedCh((subsystem, plane, bar, SiPM))
			side='left' if SiPM<8 else 'right'
			if not self.simulation:
				if fixed_ch not in self.cscintvalues:continue
				cscint=self.cscintvalues[fixed_ch]
			else: 
				cscint = MuFilter.GetConfParF(f'MuFilter/US_signalspeed_{fixed_ch}')
				if cscint==0:continue
				cscint=(cscint,0) # at the moment, no uncertainty for signal speed in simulation
			if not cscint: continue
			cscintvalues[side].append(cscint)

		average_left_cscint, average_right_cscint = sum( [ci[0] for ci in cscintvalues['left']] ) / len(cscintvalues['left']), sum( [ci[0] for ci in cscintvalues['right']] ) / len(cscintvalues['right'])
		uncertainty_sq_left, uncertainty_sq_right = sum( [ci[1]**2 for ci in cscintvalues['left']] ), sum( [ci[1]**2 for ci in cscintvalues['right']] ) 
		uncertainty_left, uncertainty_right = ROOT.TMath.Sqrt(uncertainty_sq_left), ROOT.TMath.Sqrt(uncertainty_sq_right)
	
		return (average_left_cscint, uncertainty_left), (average_right_cscint, uncertainty_right)

	def GetBarAveragesigmat(self, detID):

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

	def Makecscintdict(self, runNr, state='corrected'):
		d={}
		s=2
		for p in range(self.systemAndPlanes[s]):
			for b in range(self.systemAndBars[s]):
				for SiPM in self.systemAndSiPMs[s]:
					fixed_ch=self.MakeFixedCh((s,p,b,SiPM))
					cscint=self.Getcscint(runNr, fixed_ch=fixed_ch, state=state)
					if not cscint: continue
					else: d[fixed_ch]=cscint
		if len(d)==0: self.cscintvalues=None
		self.cscintvalues=d

	def Maketimeresolutiondict(self, runNr, state, fitmode='FWHM', histkey='dtvxpred'):

		path=f'{self.path}TimeResolution/run{runNr}/'
		res={}
		# if state=='corrected' or state=='aligned': state='corrected_tDS0-tSiPMcorrected'
		all_channels=self.GetListOfChannels(2)
		for fixed_ch in all_channels:

			filename=f'{path}timeresolution_{fixed_ch}.json'
			if not os.path.exists(filename): continue

			with open(filename, 'r') as x:
				d=json.load(x)

			if histkey not in d: 
				print(f'histkey {histkey} for {fixed_ch} not in d')
				continue 
			if state not in d[histkey]:
				print(f'state {state} for {fixed_ch} not in d[{histkey}]')
				continue
			if type(d[histkey][state])==list:
				print(f'wtf: {filename}')
			try: res[fixed_ch] = d[histkey][state][fitmode]
			except KeyError: 
				print(f'{fixed_ch} has no {state} time res for mode {fitmode}')
				continue
		self.timeresolutiondict = res

	def Writetimeresolutiondict(self, runNr, state):
		if not hasattr(self, "timeresolutiondict"): self.Maketimeresolutiondict(runNr, state)
		filename = f'{self.path}Results/run{runNr}/run{runNr}_timeresolution_{state}.json'

		# Writing JSON data to a file
		with open(filename, 'w') as json_file:
			json.dump(self.timeresolutiondict, json_file, indent=4)
		print(f'Time resolution dict written to {filename}')

	def MakeDeltatimeresolutiondict(self, runNr, state, fitmode='FWHM', histkey='dtvxpred'):
		path=f'{self.path}TimeResolution/run{runNr}/'
		self.deltatimeresolutiondict={}
		# if state=='corrected' or state=='aligned': state='corrected_tDS0-tSiPMcorrected'
		all_channels=self.GetListOfChannels(2)
		for fixed_ch in all_channels:

			filename=f'{path}timeresolutionvx_{fixed_ch}.json'
			if not os.path.exists(filename): continue

			with open(filename, 'r') as x:
				d=json.load(x)

			if histkey not in d: 
				print(f'histkey {histkey} for {fixed_ch} not in d')
				continue 
			if state not in d[histkey]:
				print(f'state {state} for {fixed_ch} not in d[{histkey}]')
				continue
			if type(d[histkey][state])==list:
				print(f'wtf: {filename}')
			try: res = d[histkey][state][fitmode]
			except KeyError: 
				print(f'{fixed_ch} has no {state} time res for mode {fitmode}')
				continue

			self.deltatimeresolutiondict[fixed_ch] = {float(k):v for k,v in res.items()}

	def GetPolyParams(self, fixed_ch, runNr='005408'):
		fname=f'{self.path}Polyparams/run{runNr}/polyparams5_{fixed_ch}.json'
		if not os.path.exists(fname): 
			return 

		with open(fname, 'r') as f:			
			alldata = json.load(f)
			data=alldata['uncorrected']

		params=[float(i) for i in data[1:13]]
		limits=[float(i) for i in data[13:15]]
		return params, limits

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

	def Gettimeresolution(self, runNr, fixed_ch, state):

		filename=f'{self.path}TimeResolution/run{runNr}/timeresolution_{fixed_ch}.json'

		if not os.path.exists(filename): return 

		with open(filename, 'r') as x:
			d=json.load(x)

		if state not in d: return 

		return d[state][0], d[state][1]		

	def GetBarsideTimeresolution(self, runNr, state, mode='FWHM'):
		
		fname = f'{self.path}Results/run{runNr}/run{runNr}_barside-timeresolutions-{state}-{mode}.json'
		with open(fname, 'r') as f:
			d = json.load(f)
		return d
	
	def CalculateBarsideTimeresolution(self, runNr, detID, side, state='corrected'):
		
		if not hasattr(self, "timeresolutiondict"): self.Maketimeresolutiondict(runNr, state)

		covar_component = self.GetCovariances(detID, side)
		SiPMresolution_component,N = self.GetResolutions(detID, side)
		
		calc_barside_resolution = np.sqrt((1/N**2)*(SiPMresolution_component + covar_component))
		return calc_barside_resolution		

	def GetResolutions(self, detID, side, state='corrected'):
		
		if side=='left':SiPMs = [0,1,3,4,6,7]
		else: SiPMs = [8,9,11,12,14,15]
		
		resolutions = []
		for SiPM in SiPMs:
			key = f'{detID}_{SiPM}'
			if not key in self.timeresolutiondict: continue
			resolution = self.timeresolutiondict[key][0]
			
			resolutions.append(resolution)
		SiPMresolution_component = sum([i**2 for i in resolutions])
		return SiPMresolution_component, len(SiPMs)

	def GetCovariances(self, detID, side):

		if not hasattr(self, "timingcovariance"): 
			self.MakeTimingCovarianceDict(self.runNr)

		if side=='left':SiPMs = [0,1,3,4,6,7]
		else: SiPMs = [8,9,11,12,14,15]
		combs = list(combinations(SiPMs, 2))
		
		covars=[]
		for i,j in combs:
			
			key = f'timingxt_{detID}_SiPMs{i}-{j}'
			if not key in self.timingcovariance: continue
			covariance = self.timingcovariance[key]
			covars.append(covariance)
		final_component = sum([2*i for i in covars])
		return final_component

	def GetSkewness(self, hist, xlow, xhigh):
		"""
		Calculate the skewness of a ROOT histogram within a specified range.
		
		Parameters:
			hist (TH1): The ROOT histogram.
			range_min (float, optional): Lower bound of the range. Use the histogram minimum if None.
			range_max (float, optional): Upper bound of the range. Use the histogram maximum if None.
		
		Returns:
			float: The skewness of the histogram in the specified range.
		"""
		# Initialize variables
		sum_weights = 0  # Total weight
		sum_x = 0        # First moment (mean numerator)
		sum_x2 = 0       # Second moment
		sum_x3 = 0       # Third moment

		# Loop over bins
		for bin_idx in range(1, hist.GetNbinsX() + 1):
			bin_center = hist.GetBinCenter(bin_idx)
			bin_content = hist.GetBinContent(bin_idx)
			
			# Only consider bins within the range
			if xlow <= bin_center <= xhigh:
				sum_weights += bin_content
				sum_x += bin_content * bin_center
				sum_x2 += bin_content * bin_center**2
				sum_x3 += bin_content * bin_center**3
		
		# Calculate mean (mu1)
		if sum_weights == 0:
			raise ValueError("No data in the specified range.")
		
		mean = sum_x / sum_weights
		
		# Calculate central moments
		mu2 = (sum_x2 / sum_weights) - mean**2
		mu3 = (sum_x3 / sum_weights) - 3 * mean * mu2 - mean**3

		# Avoid division by zero
		if mu2 <= 0:
			raise ValueError("Variance (mu2) is zero or negative; skewness is undefined.")
		
		# Calculate skewness
		skewness = mu3 / mu2**1.5
		return skewness		

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

	def Getchi2pNDF(self, runNr, fixed_ch, state, n=5):

		# iteration=0 if state=='uncorrected' else 1   
		fname=f'{self.path}chi2s/run{runNr}/chi2s{n}_{fixed_ch}.json'
		if not os.path.exists(fname): 
			print(f'no file {fname}')
			return 
		with open(fname) as jf:
			data=json.load(jf)
		chi2pNDF = data[state][0]/data[state][1]

		return chi2pNDF

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

	def GetMPV(self, runNr, fixed_ch, iteration):
		fname=f'{self.path}MPVs/run{runNr}/MPV_{fixed_ch}.csv'
		if not os.path.exists(fname): return -999.
		with open(fname, 'r') as h:
			reader=csv.reader(h)
			data=[row for row in reader]
			if data==[]: return -999. # To be investigated
			res=data[iteration-1]
		return float(res[0])

	"""
	Important note: the Analysis.correct_ToF function corrects 
	the SiPM time to x=L/2 in the physics FoR. That is not equal to the bar centre!!!! 
	"""

	def correct_ToF(self, MuFilter, fixed_ch, clock, xEx):
		detID=int(fixed_ch.split('_')[0])
		SiPM=int(fixed_ch.split('_')[-1])

		# print(f'time: {clock*6.25}, x: {xEx}')

		# Correct to the centre of the bar in the physics FoR
		MuFilter.GetPosition(detID, self.A, self.B)	
		xref = 0.5 * (self.A.x() + self.B.x())
		if not self.simulation: 
			cs = self.cscintvalues[fixed_ch]
			c_SiPM = float(cs[0])
		else: c_SiPM = MuFilter.GetConfParF(f'MuFilter/US_signalspeed_{fixed_ch}')

		# fixed_subsystem, fixed_plane, fixed_bar, fixed_SiPM = fixed
		time=clock*self.TDC2ns
		ToFcorrection=(xEx-xref)/c_SiPM
		# ToFcorrection=(xEx)/c_SiPM

		if SiPM<8: corrected_t = time + ToFcorrection
		else: corrected_t = time - ToFcorrection

		return (SiPM, corrected_t)

	def GetCorrectedTimes(self, hit, **kwargs):

		x=kwargs.get('x', 0)
		mufilter=kwargs.get('MuFilter')
		mode=kwargs.get('mode', 'aligned')

		detID=hit.GetDetectorID()

		alignedtimes=[]
		clocks, qdcs=hit.GetAllTimes(), hit.GetAllSignals()
		for i in clocks:
			SiPM, clock = i 
			fixed_ch=f'{detID}_{SiPM}'
			
			# Using qdc == -1 as a flag to only correct tof, for simulation data
			if mode == 'tof': qdc=-1
			else: qdc=self.GetChannelVal(SiPM, qdcs)

			correctedtime=self.MuFilterCorrectedTime(mufilter, fixed_ch, qdc, clock, x)
			if not correctedtime: continue
			if mode=='aligned' and not self.simulation: 
				d = self.alignmentparameters[f'{detID}_{SiPM}']
				if fixed_ch=='21000_0':
					print(f'Event number: {self.task.M.EventNumber}, reft: {self.task.reft}, twtoft: {self.task.reft-correctedtime}, d: {d[0]}')		
				correctedtime = self.task.reft - correctedtime - d[0]
			alignedtimes.append((SiPM, correctedtime))
		return alignedtimes

	def MuFilterCorrectedTime(self, MuFilter, fixed_ch, qdc, clock, x=0):
		
		time=clock*self.TDC2ns
		if not self.simulation:
			
			# To correct time: need tw params, alignment param and cscint value if correcting for signal ToF
			if not fixed_ch in self.twparameters:
				return 
			if not fixed_ch in self.alignmentparameters:
				return
			if fixed_ch not in self.cscintvalues and x!=0:
				return

			#### Correct ToF if needed.
			if x==0:
				ToFcorrectedtime = time
			else: 
				ToFcorrectedtime = self.correct_ToF(MuFilter, fixed_ch, clock, x)[1]

			# No timewalk correction if -1 passed as qdc (used for simulation data)
			if qdc==-1:
				ToFTWcorrectedtime = ToFcorrectedtime
			else:
				#### TW corrected time then ToF & TW corrected time
				twparams = self.twparameters[fixed_ch]
				twcorrection = self.correctionfunction(twparams, qdc)
				ToFTWcorrectedtime = ToFcorrectedtime+twcorrection		
			
			return ToFTWcorrectedtime

		else: 
			#### Correct ToF if needed.
			if x==0:
				ToFcorrectedtime=time
			else: 
				# if not hasattr(self.task.MuFilter.)
				if self.task.MuFilter.GetConfParF(f'MuFilter/US_signalspeed_{fixed_ch}')==0:
					return
				ToFcorrectedtime=self.correct_ToF(fixed_ch, clock, x)[1]
			return ToFcorrectedtime

	def GetTimingDiscriminant(self, hits):

		US1hits=[h.GetDetectorID() for h in hits if all([h.GetDetectorID()//10000==2, self.parseDetID(h.GetDetectorID())[1]==0])]

		# If no US1 hits found, disregard event
		if len(US1hits)==0: return -999

		# If 1 US1 hit found, ideal
		elif len(US1hits)==1:
			x=US1hits[0]
			tmp={h.GetDetectorID():h for h in hits}
			us1hit=tmp[x]

		# If multiple US1 hits found, take hit with lowest DOCA to extrapolated DS tracklet
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

		# averageUS1time=self.GetAverageTime(us1hit, correctTW=True)
		# if not averageUS1time: return -998 

		DS3Haverage=self.GetDSHaverage(hits, mode='timingdiscriminant')

		# Manual determination of timing disc.
		res=[]
		us1correctedtimes = self.GetCorrectedTimes(us1hit)
		detID = us1hit.GetDetectorID()
		for x in us1correctedtimes:
			SiPM, ctime = x 
			fixed_ch = f'{detID}_{SiPM}'
			if not fixed_ch in self.alignmentparameters: continue
			d = self.alignmentparameters[fixed_ch]
			res.append(DS3Haverage - (ctime - d[0]))
		if len(res)==0: return -999
		return sum(res)/len(res)

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
			elif self.timealignment=='new': filename=f'{self.path}rootfiles/run005408/SelectionCriteria.root'
			elif self.timealignment=='new+LHCsynch': filename=f'{self.path}rootfiles/run005999/SelectionCriteria.root'

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
		fname=f'{self.path}Alignmentparams/run{runNr}/alignmentparams_{self.refsysname}_{fixed_ch}.csv'
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

		if state=='aligned': histname=f'dt_{fixed_ch}_{self.refsysname}aligned'
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

	def MakeTWCorrectionDict(self, withErrors=False):

		d={}
		for s in (2,): # Only make dict for US
			for p in range(self.systemAndPlanes[s]):
				for b in range(self.systemAndBars[s]):
					for SiPM in self.systemAndSiPMs[s]:
						fixed_ch=self.MakeFixedCh((s,p,b,SiPM))
						tmp=self.GetPolyParams(fixed_ch)
						if not tmp: 
							print(f'No tw params for {fixed_ch}')
							continue
						else: paramsAndErrors=tmp[0]
						if not withErrors: params=self.correctionparams(paramsAndErrors)
						else: params=paramsAndErrors
						d[fixed_ch]=params
		if len(d)==0: self.twparameters=None
		self.twparameters=d

	def WriteTWParamDict(self):
		if not hasattr(self, "twparameters"): self.MakeTWCorrectionDict(self.timealignment)

		twparamsdir = f'{self.path}/Polyparams/run{self.runNr}/'
		twparamsfilename=twparamsdir+f'twparams.json'	
		with open(twparamsfilename, 'w') as jf: 
			json.dump(self.twparameters, jf)
		print(f'TW param dict written to {twparamsfilename}')

	### Make dictionary of the alignment parameter determined as the truncated y-mean of tw-corr (tds0 - tSiPM)
	def MakeAlignmentParameterDict(self):

		# Update stored string for alignment params stored
		# self.timealignment=alignment

		# if alignment=='old': run=str(5097).zfill(6)
		# elif alignment=='new': run=str(5408).zfill(6)
		# elif alignment=='new+LHCsynch': run=str(5999).zfill(6)

		d={}
		for s in (2,): # Just US
			for p in range(self.systemAndPlanes[s]):
				for b in range(self.systemAndBars[s]):
					for SiPM in self.systemAndSiPMs[s]:
						fixed_ch=self.MakeFixedCh((s,p,b,SiPM))
						correction=self.GetAlignmentParameters(self.runNr, fixed_ch)
						if not correction: continue
						d[fixed_ch]=correction
		if len(d)==0: self.alignmentparameters=None
		self.alignmentparameters=d

	def WriteAlignmentParamDict(self):
		if not hasattr(self, "alignmentparameters"): self.MakeAlignmentParameterDict()

		alignmentparamsdir = f'{self.path}/Alignmentparams/run{self.runNr}/'
		alignmentparamsfilename=alignmentparamsdir+f'alignmentparams_{self.refsysname}.json'	
		with open(alignmentparamsfilename, 'w') as jf: 
			json.dump(self.alignmentparameters, jf)
		print(f'Alignment param dict written to {alignmentparamsfilename}')

	def MakeTimingCovarianceDict(self, runNr, mode='truncated'):

		if mode=='truncated': filename=f'{self.path}TimingCovariance/run{runNr}/run{runNr}_truncatedcovariance.json'
		elif mode=='full': filename=f'{self.path}TimingCovariance/run{runNr}/timingcovariance.json'
		else: return

		with open(filename, 'r') as x:
			self.timingcovariance = json.load(x)

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
		mpvfilename=mpvpath+f'MPVs.json'
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

	def MakePlaneWiseBarMapThing(self, MuFilter):

		d={}
		for plane in range(5):
			for bar in range(10):
				detID = int(f'2{plane}00{bar}')
				MuFilter.GetPosition(detID, self.A, self.B)
				d[detID] = [self.A.x(), self.B.x(), self.A.y(), self.B.y()]
		for plane in range(4):
			if plane!=3:bars=range(120)
			else: bars=range(60,120)
			for bar in bars:
				bar=str(bar).zfill(3)
				detID=int(f'3{plane}{bar}')
				MuFilter.GetPosition(detID, self.A, self.B)
				d[detID] = [self.A.x(), self.B.x(), self.A.y(), self.B.y()]				

		for f in (f'{self.path}BarPositions.json', f'/eos/user/a/aconsnd/SWAN_projects/numuInvestigation/data/BarPositions.json', f'/eos/user/a/aconsnd/SWAN_projects/Simulation/data/BarPositions.json'):
			with open(f, 'w') as jf:
				json.dump(d, jf)
		print(f'Bar positions written to three locations')
		return d
