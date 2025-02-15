#!/usr/bin/env python3

import os, json
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

# Set the default font size for titles
plt.rcParams['figure.titlesize'] = 18
# Set the default font size for titles
plt.rcParams['axes.titlesize'] = 16
# Set the default font size for labels
plt.rcParams['axes.labelsize'] = 16
fontsize=16

class TimeWalk(object):
	def __init__(self, M):

		self.options = M.options
		self.M = M 
		self.intensity_colour='C1'
		self.qdc_colour='C0'
		self.colours = ['blue', 'green', 'red', 'cyan', 'magenta', 'yellow', 'tab:orange', 'tab:purple']
		self.GetQDCvTimestamp()
		self.dts={'TI18':{}, 'laser':{}}

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
		
		d=f'{self.M.sndswpath}analysis-plots/timewalk/'
		os.makedirs(d, exist_ok=True)
		
		plotlocation = d+f'SiPM{self.M.SiPM}_{self.M.titledict[self.M.mode]}.png'
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

		for i in range(2): axes[i].legend(fontsize=fontsize)

		d=f'{self.M.sndswpath}analysis-plots/timewalk/'
		os.makedirs(d, exist_ok=True)
		
		plotlocation = d+f'SiPM{self.M.SiPM}_{self.M.titledict[self.M.mode]}.png'
		fig.savefig(plotlocation, bbox_inches='tight')
		print(f'Plot saved to {plotlocation}')
	
	def CompareRunTypes(self, style='overlay'):
		if style=='overlay':
			fig, ax = plt.subplots(figsize=(10,6))
			ax.tick_params(axis='both', which='major', labelsize=fontsize)

			ax.set_title(f"Timestamp v QDC for SiPM {self.M.SiPM}, bar {self.M.bar}\nwith different illumination modes", fontsize=fontsize)
			ax.set_ylabel('Timestamp [ns]', fontsize=fontsize)
			ax.set_xlabel('QDC [a.u]', fontsize=fontsize)
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
			ax.legend(loc=(0.5, 0.2), fontsize=fontsize)
				
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

		d=f'{self.M.sndswpath}analysis-plots/timewalk/illuminations/'
		os.makedirs(d, exist_ok=True)
		
		plotlocation = d+f'compare-run-types_SiPM{self.M.SiPM}.png'
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

		# Invert trend just for better comparison with laser, done for DPG2024
		tmp = [(twfit1Dhist.GetBinLowEdge(i), -twfit1Dhist.GetBinContent(i), twfit1Dhist.GetBinError(i)) for i in range(twfit1Dhist.GetNbinsX())] 
		x = [i for i in tmp if i[1]!=0]

		TI18_qdcs, TI18_times, TI18_errors = zip(*x)

		return TI18_qdcs, TI18_times, TI18_errors, fixed_ch

	def TI18_laser_comparison(self):

		# I will add loops for comparable SiPMs
		idx=0

		if len(self.data)==0: return
		qdcs, times = zip(*self.data.values())

		lasermip_idx = next((i for i,x in enumerate(qdcs) if x >= 5), None) # pass generator object to next, return None if all < 5 which will never happen
		lasermip = qdcs[lasermip_idx], times[lasermip_idx] 
		laser100_idx = min(range(len(qdcs)), key=lambda i:abs(qdcs[i]-100)) # return value closest to 100
		laser100 = qdcs[laser100_idx], times[laser100_idx]

		TI18_qdcs, TI18_times, TI18_errors, fixed_ch = self.GetTI18Timewalk(
			plane=self.options.tw_comparison_plane,
			side=self.options.tw_comparison_side
			)
		TI18mip_idx = next((i for i,x in enumerate(TI18_qdcs) if x==lasermip[0]), None)
		TI18mip = TI18_qdcs[TI18mip_idx], TI18_times[TI18mip_idx]
		TI18100_idx = next((i for i,x in enumerate(TI18_qdcs) if x==laser100[0]), None)
		TI18100 = TI18_qdcs[TI18100_idx], TI18_times[TI18100_idx]	

		if any([TI18100==None, TI18mip==None]): return
		self.dts['TI18'][self.M.SiPM]=abs(TI18100[1]-TI18mip[1])
		self.dts['laser'][self.M.SiPM]=abs(laser100[1]-lasermip[1])

		laserlabel = f't({laser100[0]}) - t({lasermip[0]}) = {round(laser100[1]-lasermip[1], 2)} ns'
		TI18label = f't({TI18100[0]}) - t({TI18mip[0]}) = {round(TI18100[1]-TI18mip[1], 2)} ns'

		# Get graph points from timewalk_{fixed_ch}.root
		fig, axes = plt.subplots(1,2, figsize=(14,6))

		fig_title=fig.suptitle(f'Comparison of timewalk observed in TI18 with laser data SiPM {self.M.SiPM}\n{self.M.titledict[self.M.mode]} illumination')	
		for i in range(2): 
			axes[i].tick_params(axis='both', which='major', labelsize=fontsize)
			axes[i].set_xlabel('QDC [a.u]', position=(1,0), horizontalalignment='right', fontsize=fontsize)
			axes[i].grid(which='major', axis='both', linestyle='--')
			
		axes[0].set_ylabel('$-(t_{0}^{DS} - t_{SiPM})$ [ns]', position=(0,1), horizontalalignment='right', fontsize=fontsize)
		axes[0].set_title(f'TI18: {self.M.LaserAna.MakeHumanReadableFixedCh(fixed_ch)}', fontsize=fontsize)
		axes[1].set_ylabel('time [ns]', position=(0,1), horizontalalignment='right', fontsize=fontsize)
		axes[1].set_title('Laser', fontsize=fontsize)

		# TI18 data
		axes[0].errorbar(TI18_qdcs, [i for i in TI18_times], color="C1", 
			yerr=TI18_errors, fmt='o', label=TI18label) 
		axes[0].legend(fontsize=fontsize)

		# Laser data
		axes[1].errorbar(qdcs, times, color=self.colours[idx], 
			xerr=[d[0] for d in self.errors.values()], yerr=[d[1] for d in self.errors.values()], 
			label=laserlabel)
		axes[1].legend(fontsize=fontsize)

		points_to_label={0:{
							TI18mip_idx:[TI18_qdcs[TI18mip_idx], TI18_times[TI18mip_idx]],TI18100_idx:[TI18_qdcs[TI18100_idx], TI18_times[TI18100_idx]]
						}, 
						1:{
							lasermip_idx:[qdcs[lasermip_idx], times[lasermip_idx]],
							laser100_idx:[qdcs[laser100_idx], times[laser100_idx]]
						}
		}

		for idx,ax in enumerate(axes):
			for pt, d in points_to_label[idx].items():
				ax.annotate(f'({round(d[0],1)}, {round(d[1],1)})', (d[0], d[1]), va='bottom', ha='left',
				# xytext=(0.5,0.5),
				fontsize=fontsize)

		# Make directory if not already there
		d=f'{self.M.sndswpath}analysis-plots/timewalk/TI18-comparison/'
		os.makedirs(d, exist_ok=True)

		plt.tight_layout()
		plotlocation = d+f'SiPM_{self.M.SiPM}_{self.M.titledict[self.M.mode]}.png'
		fig.savefig(plotlocation, bbox_inches='tight')
		print(f'Plot saved to {plotlocation}')

	def LoadTI18comparisondata(self):
		
		d = f'{self.M.sndswpath}analysis-plots/timewalk/TI18-comparison/'
		picklelocation= d+f'{self.M.titledict[self.M.mode]}-comparison.json'
		with open(picklelocation, 'r') as f:
			self.dt = json.loads(f.read())

	def PlotSystemTI18LaserDifference(self):

		# fig, axes = plt.subplots(2,1, figsize=(8, 12))
		fig, axes = plt.subplots(2, 1, sharex=True, figsize=(8, 8))
		fig_title=fig.suptitle(f'Comparison of time walk between laser and TI18 for all SiPMs\n{self.M.titledict[self.M.mode]} illumination')

		for i in axes.flatten():
			
			i.grid(which='both', axis='y', linestyle='--')
			[i.axvline(g*8, linestyle='--', color='grey') for g in range(1,10)]
			i.tick_params(axis='both', which='major', labelsize=fontsize)
		axes[1].set_xlabel('SiPM number', position=(1,0), horizontalalignment='right', fontsize=fontsize)
		axes[0].set_ylabel('|dt| [ns]', position=(0,1), horizontalalignment='right', fontsize=fontsize)
		axes[0].set_title(f'Degree of timewalk observed over fixed QDC range', fontsize=fontsize)
		axes[1].set_ylabel('Timewalk laser / TI18', position=(0,1), horizontalalignment='right', fontsize=fontsize)
		# axes[1].set_title('Ratio of timewalk observed between TI18 and laser', fontsize=fontsize)

		for i, x in enumerate((('laser', 'o'), ('TI18', 'D'))):
			axes[0].errorbar([int(i) for i in self.dts[x[0]].keys()],self.dts[x[0]].values(), label=f'{x[0]} data', marker=x[1])
		axes[0].legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon=True, fontsize=fontsize)

		plane = self.options.tw_comparison_plane
		side = self.options.tw_comparison_side

		d_data = {int(SiPM):abs(self.dts['TI18'][SiPM]/self.dts['laser'][SiPM]) for SiPM in self.dts['laser'].keys()}
		axes[1].errorbar(d_data.keys(), d_data.values())

		axes[1].text(
				0.65, 0.2, 
				f"TI18 data from plane {plane+1}, {side} side", 
				color='black', fontsize=fontsize, ha='left', va='center', transform=axes[1].transAxes,
    			bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.5', alpha=0.5)
				)


		d = f'{self.M.sndswpath}analysis-plots/timewalk/TI18-comparison/'

		plotlocation = d+f'SystemTimewalkComparison_{self.M.titledict[self.M.mode]}-plane{plane}side{side}.png'
		plt.tight_layout()
		fig.savefig(plotlocation, bbox_inches='tight')
		print(f'Timewalk comparison plot saved to {plotlocation}')

		jsonlocation = d+f'{self.M.titledict[self.M.mode]}-comparison.json'
		with open(jsonlocation, 'w') as f:
			json.dump(self.dts, f)
		print(f'Timewalk comparison data written to {jsonlocation}')
