#!/usr/bin/env python3

import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

# Set the default font size for titles
plt.rcParams['figure.titlesize'] = 18
# Set the default font size for titles
plt.rcParams['axes.titlesize'] = 16
# Set the default font size for labels
plt.rcParams['axes.labelsize'] = 16

class TimeWalk(object):
	def __init__(self, M):

		self.options = M.options
		self.M = M 
		self.intensity_colour='C1'
		self.qdc_colour='C0'
		self.colours = ['blue', 'green', 'red', 'cyan', 'magenta', 'yellow', 'tab:orange', 'tab:purple']
		self.GetQDCvTimestamp()

	def GetQDCvTimestamp(self):

		timewalk_runs = self.M.runs[self.M.mode]
		self.data={}
		self.errors={}
		self.intensities=[]

		for intensity in timewalk_runs:
			runNr = timewalk_runs[intensity]
			timestamp = self.M.GetHistogramValue(runNr, 'timestamp')
			if not timestamp: continue
			qdc = self.M.GetHistogramValue(runNr, 'QDC')

			if any([not timestamp, not qdc]): 
				# print(f'Run {runNr} has no SiPM {self.M.SiPM} plots made')
				continue
			
			timing_ns = [i*self.M.TDC2ns for i in timestamp]
			self.data[runNr]=(qdc[0],timing_ns[0])
			self.errors[runNr]=(qdc[1],timing_ns[1])
			self.intensities.append(intensity)
		# return intensities, data, errors

	def SetMode(self, mode):
		self.M.SetMode(mode)
		self.GetQDCvTimestamp()

	def SetSiPM(self, SiPM):
		self.M.SetSiPM(SiPM)
		self.GetQDCvTimestamp()

	def TimeWalk(self):

		if len(self.data)==0: return
		qdcs, timestamps = zip(*self.data.values())

		labels=[]
		for i in range(len(self.data.keys())):
			if i==0:labels.append(f'run {list(self.data.keys())[i]}')
			else:
				dt = self.data[list(self.data.keys())[i-1]][1]-self.data[list(self.data.keys())[i]][1]
				dqdc = self.data[list(self.data.keys())[i-1]][0]-self.data[list(self.data.keys())[i]][0]
				if dqdc==0: continue
				label=f'run {list(self.data.keys())[i]}'
				labels.append(label)

		fig, qdc_ax = plt.subplots(figsize=(10,8))

		qdc_ax.errorbar(qdcs, timestamps, 
			xerr=[d[0] for d in self.errors.values()], yerr=[d[1] for d in self.errors.values()], 
			label=f'SiPM {self.M.SiPM}')		
		intensity_ax=qdc_ax.twiny()
		intensity_ax.set_label('Intensity values')
		intensity_ax.plot(self.intensities, timestamps, color=self.intensity_colour,marker='+',linestyle='--', label='Intensity values')
		qdc_ax.tick_params(labelright=True)
		fig.tight_layout()

		for label, i,j in zip(labels, qdcs, timestamps):    
			qdc_ax.annotate(label, (i, j), textcoords="offset points", xytext=(0,10), ha='center')

		plt.minorticks_on()
		plt.gca().xaxis.set_minor_locator(AutoMinorLocator(n=5))  
		plt.gca().yaxis.set_minor_locator(AutoMinorLocator(n=5))
        
		qdc_ax.set_xlabel('QDC [a.u]')
		intensity_ax.set_xlabel('Laser intensity [%]')
		qdc_ax.set_ylabel('Time stamp [ns]')

		# Combine legend entries for both axes
		lines, labels = qdc_ax.get_legend_handles_labels()
		lines2, labels2 = intensity_ax.get_legend_handles_labels()
		qdc_ax.legend(lines + lines2, labels + labels2, loc='upper right')        

		qdc_ax.grid(True, linestyle='--', alpha=0.7, color='gray')
		qdc_ax.set_title(f'Time walk measured with laser for SiPM {self.M.SiPM}\n{self.M.titledict[self.M.mode]} illumination')
		plotlocation = f'{self.M.sndswpath}analysis-plots/timewalk/SiPM{self.M.SiPM}_{self.M.titledict[self.M.mode]}.png'
		plt.savefig(f'{plotlocation}', bbox_inches='tight')
		print(f'SiPM {self.M.SiPM} timewalk figure saved to {plotlocation}')

	def OverlaySiPMs(self, bar):
		sipmtypekey={0:'large', 1:'small'}
		SiPMs = list(range((bar-1)*8+1, bar*8+1))

		fig, axes = plt.subplots(1,2, figsize=(16,8))

		fig_title=fig.suptitle(f'Overlay of large and small SiPM timewalk plots\n{self.M.titledict[self.M.mode]} illumination')	
		# fig_title.set_y(0.1)
		for i in range(2):
			axes[i].set_title(f"Timestamp v QDC for {sipmtypekey[i]} SiPMs in bar {bar}")
			axes[i].set_xlabel('QDC [a.u]')
			axes[i].set_ylabel('Timestamp [ns]')

		for idx,SiPM in enumerate(SiPMs):	

			# intensities, data, errors = self.GetQDCvTimestamp()
			self.SetSiPM(SiPM)
			if len(self.intensities)==0: continue

			qdcs, timestamps = zip(*self.data.values())

			if self.M.SiPMsize=='large': 
				axes[0].errorbar(qdcs, timestamps, color=self.colours[idx], 
					xerr=[d[0] for d in self.errors.values()], yerr=[d[1] for d in self.errors.values()], 
					label=f'SiPM {self.M.SiPM}')			

			if self.M.SiPMsize=='small': 
				axes[1].errorbar(qdcs, timestamps, color=self.colours[idx], 
					xerr=[d[0] for d in self.errors.values()], yerr=[d[1] for d in self.errors.values()], 
					label=f'SiPM {self.M.SiPM}')				

		for i in range(2): axes[i].legend()

		plotlocation = f'{self.M.sndswpath}analysis-plots/timewalk/bar{bar}_{self.M.titledict[self.M.mode]}.png'
		fig.savefig(plotlocation, bbox_inches='tight')
		print(f'Plot saved to {plotlocation}')
	
	def CompareRunTypes(self, style='overlay'):
		if style=='overlay':
			fig, ax = plt.subplots()

			ax.set_title(f"Timestamp v QDC for SiPM {self.M.SiPM}, bar {self.M.bar}\nwith different illumination modes")
			ax.set_ylabel('Timestamp [ns]')
			ax.set_xlabel('QDC [a.u]')
			ax.grid(which='major',axis='both',linestyle='--')
			for idx,runtype in enumerate(self.M.runs.keys()):
				
    			# For SiPMs with no 3-bar illumination data:
				if len(self.M.runs[runtype])==0: continue

				self.M.mode=runtype
				self.GetQDCvTimestamp()
    
				if len(self.intensities)==0: continue
				qdcs, timestamps = zip(*self.data.values())    
    
				if runtype=='SiPM':	label = f"Single SiPM {self.M.SiPM}"
				elif runtype==1:	label = f"1-bar illuminated"
				else:				label = f"{runtype}-bars illuminated"    

				if self.M.SiPMsize=='large': 
					ax.errorbar(qdcs, timestamps, color=self.colours[idx], 
						xerr=[d[0] for d in self.errors.values()], yerr=[d[1] for d in self.errors.values()], 
						label=label)			

				if self.M.SiPMsize=='small': 
					ax.errorbar(qdcs, timestamps, color=self.colours[idx], 
						xerr=[d[0] for d in self.errors.values()], yerr=[d[1] for d in self.errors.values()], 
						label=label)
			ax.legend()
				
		if style=='subplots':
			fig, axes = plt.subplots(2,2, figsize = (12, 12))
			fig_title=fig.suptitle(f"QDC v timestamps for SiPM {self.M.SiPM}\nwith different illumination modes")

			for idx,runtype in enumerate(self.M.runs.keys()):
				self.M.mode=runtype
				# current, qdcs, errors = self.GetQDCvIntensity()
				if runtype=='SiPM':	axes[0,0].set_title(f"Single SiPM")
				elif runtype==1:	axes[0,1].set_title(f"1-bar illuminated")
				else:				axes[idx//2, idx%2].set_title(f"{runtype}-bars illuminated")

				axes[idx//2, idx%2].set_xlabel('QDC [a.u]')
				axes[idx//2, idx%2].set_ylabel('Timestamp [ns]')
				axes[idx//2, idx%2].grid(which='major', axis='both', linestyle='--')
				axes[idx//2, idx%2].errorbar(qdcs, timestamps, color=self.colours[idx], 
						xerr=[d[0] for d in self.errors.values()], yerr=[d[1] for d in self.errors.values()], 
						label=f'SiPM {self.M.SiPM}')

		plotlocation = f'{self.M.sndswpath}analysis-plots/timewalk/illuminations/{style}-SiPM{self.M.SiPM}.png'
		fig.savefig(plotlocation, bbox_inches='tight')
		print(f'Plot saved to {plotlocation}')		

     
	def GetTI18Timewalk(self, plane=0, side='left'):

		# Need ROOT to get TI18 data
		from ROOT import TFile, TPad, TH1D, gROOT

		bar, fSiPM = self.M.LaserAna.SiPM2BarAndPosition(self.M.SiPM)
		if side=='right': fSiPM = fSiPM+8
		fSiPM -= 1
		bar-=1

		fixed_ch = self.M.LaserAna.MakeFixedCh((2,plane,bar,fSiPM))
		# Put TI18 TW path in self.M
		filename = f'/afs/cern.ch/work/a/aconsnd/Timing-physics2022/rootfiles/run005408/timewalk_{fixed_ch}.root'
		f=TFile.Open(filename, 'read')
		twfitcanvas=f.Get(f'twfitcanvas_{fixed_ch}')
		twfitpad = twfitcanvas.GetPrimitive(f"twfitpad_{fixed_ch}")
		twfit1Dhist = twfitpad.GetPrimitive(f"dtvqdc_{fixed_ch}_uncorrected_px").Clone()
		twfit1Dhist.SetDirectory(gROOT)
		f.Close()

		# self.twhist=twfit1Dhist

		tmp = [(twfit1Dhist.GetBinCenter(i), twfit1Dhist.GetBinContent(i), twfit1Dhist.GetBinError(i)) for i in range(twfit1Dhist.GetNbinsX())]
		x = [i for i in tmp if i[1]!=0]

		TI18_qdcs, TI18_timestamps, TI18_errors = zip(*x)

		return TI18_qdcs, TI18_timestamps, TI18_errors, fixed_ch

	def TI18_laser_comparison(self):

		# I will add loops for comparable SiPMs
		idx=0

		if len(self.data)==0: return
		qdcs, timestamps = zip(*self.data.values())

		# Find delta t between first qdc>5 and last qdc<100
		lasermip, laser100 = None,None
		for i in range(len(qdcs)):
			if qdcs[i]>=5 and lasermip==None:
				lasermip = (qdcs[i], timestamps[i])
			if qdcs[i]>=100:
				laser100 = (qdcs[i-1], timestamps[i-1])
				break 				
		
		TI18_qdcs, TI18_timestamps, TI18_errors, fixed_ch = self.GetTI18Timewalk()

		# Find delta t between first qdc>5 and last qdc<100
		TI18mip, TI18100 = None,None
		for i in range(len(TI18_qdcs)):
			if TI18_qdcs[i]>=5 and TI18mip==None:
				TI18mip = (TI18_qdcs[i], TI18_timestamps[i])
			if TI18_qdcs[i]>=100:
				TI18100 = (TI18_qdcs[i-1], TI18_timestamps[i-1])
				break 		

		laserlabel = f't({laser100[0]}) - t({lasermip[0]}) = {round(laser100[1]-lasermip[1], 2)} ns'
		TI18label = f't({TI18100[0]}) - t({TI18mip[0]}) = {round(TI18100[1]-TI18mip[1], 2)} ns'

		# Get graph points from timewalk_{fixed_ch}.root
		fig, axes = plt.subplots(1,2, figsize=(14,8))

		fig_title=fig.suptitle(f'Comparison of timewalk observed in TI18 with laser data SiPM {self.M.SiPM}\n{self.M.titledict[self.M.mode]} illumination')	
		for i in range(2): 
			axes[i].set_xlabel('QDC [a.u]', position=(1,0), horizontalalignment='right')
			axes[i].grid(which='major', axis='both', linestyle='--')
		axes[0].set_ylabel('$t_{0}^{DS} - t_{SiPM} [ns]$', position=(0,1), horizontalalignment='right')
		axes[0].set_title(f'TI18: {self.M.LaserAna.MakeHumanReadableFixedCh(fixed_ch)}')
		axes[1].set_ylabel('timestamp [ns]', position=(0,1), horizontalalignment='right')
		axes[1].set_title('Laser')

		# TI18 data
		axes[0].errorbar(TI18_qdcs, TI18_timestamps, color="C1",
			yerr=TI18_errors, fmt='o', label=TI18label)
		axes[0].legend()

		# Laser data
		axes[1].errorbar(qdcs, timestamps, color=self.colours[idx], 
			xerr=[d[0] for d in self.errors.values()], yerr=[d[1] for d in self.errors.values()], 
			label=laserlabel)
		axes[1].legend()

		plotlocation = f'{self.M.sndswpath}/analysis-plots/timewalk/TI18-comparison/SiPM_{self.M.SiPM}_{self.M.titledict[self.M.mode]}.png'
		fig.savefig(plotlocation, bbox_inches='tight')
		print(f'Plot saved to {plotlocation}')