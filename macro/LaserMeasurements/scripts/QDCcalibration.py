#!/usr/bin/env python3

import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import numpy as np
import json
from scipy.optimize import curve_fit
import csv, os
plt.rcParams['font.size'] = 16
plt.rcParams['lines.markersize'] = 8

class QDCcalibration(object):
	def __init__(self, M):

		self.options = M.options
		self.M = M
		self.laser_mode = self.options.laser_mode

		self.SiPMconfigs=self.GetSiPMconfigs()
		if self.laser_mode.find('ngamma')==0: self.GetPhotonConversionConstants()
		self.lasermode_displaydict = {'uncalibrated':'Uncalibrated laser intensity [%]', 'calibrated':'Calibrated laser intensity [%]', 
            						'adc':'Average adc [a.u]', 'ngamma_incident': 'Incident photons, $n_{\gamma}$',
									'ngamma_detected': 'Detected photons, $n_{\gamma}$', 'Npe':"N_{pe}"} 
		self.int_title = self.lasermode_displaydict[self.laser_mode]		
		self.colours = ['blue', 'green', 'red', 'cyan', 'magenta', 'yellow', 'tab:orange', 'tab:purple']
		self.calibrated_intensities=self.M.calibrated_intensities
		self.saturation_params={}
		self.saturation_covariance = {}
		self.linear_params={}
		self.linear_covariance = {}
		self.GetQDCvIntensity()
		
	def GetPhotonConversionConstants(self):
		self.area, self.PDE, self.pixels, self.gain = self.SiPMconfigs[self.M.SiPMsize].values()

	def GetSiPMconfigs(self):
		configfile = f"{self.M.eospath}/configuration/SiPMparameters.json"
		with open(configfile, 'r') as f:
			d = json.load(f) 
		return d		

	def saturationfunction(self, x, a, b, c, d):
		return a * x - b * np.exp(-c * x) + d
	def linearfunction(self, x, a, b):	
		return a*x + b
		
	def SetSiPM(self, SiPM):
		self.M.SetSiPM(SiPM)
		self.GetQDCvIntensity()

	def GetQDCvIntensity(self):
		mpv_runs = self.M.runs[self.M.mode]
		data={}
		intensities=[] 

		for intensity in mpv_runs:
			runNr = mpv_runs[intensity]
			qdc = self.M.GetHistogramValue(runNr, 'QDC')
			if not qdc: continue

			data[runNr]=qdc
			intensities.append(intensity)
        
		labels=[]
		for i in range(len(data.keys())):
			labels.append(f'run {list(data.keys())[i]}')

		self.qdcs=[i[0] for i in list(data.values())]
		self.errors=[abs(i[1]) for i in list(data.values())]

		for i in range(len(self.errors)): 
			if self.errors[i]<0: 
				print(f'Run {mpv_runs[intensities[i]]} has a negative dQDC')

		if self.laser_mode.find('ngamma')==-1 or self.laser_mode != 'N_pe': df_key = self.lasermode_displaydict[self.laser_mode]
		else: df_key='Incident photons / mm_sq'

		self.current = self.calibrated_intensities[self.calibrated_intensities["Laser intensity [%]"].isin(intensities)][df_key]
		if self.laser_mode.find('ngamma')!= -1: self.current = self.current*self.area
		if self.laser_mode.find('detected')!= -1: self.current = self.current*self.PDE
		if self.laser_mode == 'Npe': self.current = self.current*self.area*self.PDE*self.gain
	def MPVprogression(self):

		if len(self.current)==0:
			print(f"No data for SiPM {self.M.SiPM}")
			return

		fig, ax = plt.subplots(figsize=(10,6))
		ax.set_xlim(0,max(self.current))
		ax.errorbar(self.current.array, self.qdcs, yerr=self.errors, fmt='o', linestyle='-', markersize=8, label=f'SiPM {self.M.SiPM}')

		if self.laser_mode.find('gamma')!=-1:
			fracpixels_ax=ax.twiny()
			fracpixels_ax.set_label('$n_{\gamma} / n_{pixels}$')
			# fracpixels_ax.plot(self.current.array/self.pixels, self.qdcs, colour='C1' label='$n_{\gamma} / n_{pixels}$')  
			self.frac = 1-np.exp(-self.current.array * self.PDE * 1/self.pixels)
			fracpixels_ax.plot(self.frac, self.qdcs, alpha=0)  
			fracpixels_ax.set_xlabel('$n_{\gamma} / n_{pixels}$')

			# Combine legend entries for both axes
			lines, labels = ax.get_legend_handles_labels()
			lines2, labels2 = fracpixels_ax.get_legend_handles_labels()

		elif self.laser_mode=='adc':
			self.DetermineQDCoffset()
			saturation_params = self.saturation_params[self.M.SiPM]
			linear_params = self.linear_params[self.M.SiPM]

			N=np.linspace(0, max(self.current), 1000)
			y_saturation_fit = self.saturationfunction(N, *saturation_params)
			y_linear_fit = self.linearfunction(N, *linear_params)

			ax.plot(N, np.array(y_saturation_fit), label=f'ax - b(exp^(-cx)) + d ', color='orange')
			ax.plot(N, np.array(y_linear_fit), label=f'ax + b (fitted to first 6 data points)', color='green')	
			if self.M.SiPMsize=='large':text_position = (500, 50)
			else: text_position = (1000, -2.5)

			if self.options.fitmode=='linear': ax.text(*text_position, f'linear offset = {round(linear_params[1], 1)} QDC', fontsize=14)
			elif self.options.fitmode=='inv_exp': ax.text(*text_position, f'inv. exp. offset = {round(saturation_params[-1], 1)} QDC', fontsize=14)

		elif self.laser_mode=='Npe':
			self.DetermineQDCoffset()
			saturation_params = self.saturation_params[self.M.SiPM]
			linear_params = self.linear_params[self.M.SiPM]

			N=np.linspace(0, max(self.current), 1000)
			y_saturation_fit = self.saturationfunction(N, *saturation_params)
			y_linear_fit = self.linearfunction(N, *linear_params)

			ax.plot(N, np.array(y_saturation_fit), label=f'ax - b(exp^(-cx)) + d ', color='orange')
			ax.plot(N, np.array(y_linear_fit), label=f'ax + b (fitted to first 6 data points)', color='green')	
			if self.M.SiPMsize=='large':text_position = (500, 50)
			else: text_position = (1000, -2.5)

			if self.options.fitmode=='linear': ax.text(*text_position, f'linear offset = {round(linear_params[1], 1)} QDC', fontsize=14)
			elif self.options.fitmode=='inv_exp': ax.text(*text_position, f'inv. exp. offset = {round(saturation_params[-1], 1)} QDC', fontsize=14)			
  
		plt.minorticks_on()
		plt.gca().xaxis.set_minor_locator(AutoMinorLocator(n=5))  
		plt.gca().yaxis.set_minor_locator(AutoMinorLocator(n=5))

		ax.set_ylabel('QDC [a.u]')
		
		ax.set_xlabel(self.int_title)

		if self.laser_mode.find('gamma')!=-1: ax.legend(lines + lines2, labels + labels2, loc='lower right')
		else: ax.legend()

		ax.grid(True, linestyle='--', alpha=0.7, color='gray')
		ax.set_title(f'SiPM {self.M.SiPM} QDC measured as a function of {self.int_title.lower()}\n{self.M.titledict[self.M.mode]} illumination')
		
		# Create directories if not existing already
		d=f'{self.M.sndswpath}/analysis-plots/qdc_calibration/'
		os.makedirs(d, exist_ok=True)

		plotlocation = d+f'SiPM{self.M.SiPM}_{self.M.titledict[self.M.mode]}_{self.laser_mode}.png'
		plt.savefig(f'{plotlocation}', bbox_inches='tight')
		print(f'SiPM {self.M.SiPM} qdc v intensity figure saved to {plotlocation}')

	def DetermineQDCoffset(self):

		if len(self.current)==0:
			print(f'No data for SiPM {self.M.SiPM}')
			return
		# quadratic_fit = np.polyfit(self.current, self.qdcs, 2)
		# linear_component, offset = quadratic_fit[1], quadratic_fit[2]
		p0s = [1, 1, 2, 0]
		lbs, ubs = [0,0,0,-25], [5,5,5,25]
		bounds = [lbs, ubs]

		saturation_params, saturation_covariance = curve_fit(self.saturationfunction, self.current.array, self.qdcs, p0=p0s, bounds=bounds)
		linear_params, linear_covariance = curve_fit(self.linearfunction, self.current.array[:6], self.qdcs[:6])

		self.saturation_params[self.M.SiPM] = saturation_params 
		self.saturation_covariance[self.M.SiPM] = saturation_covariance

		self.linear_params[self.M.SiPM] = linear_params 
		self.linear_covariance[self.M.SiPM] = linear_covariance  
		# y_fit = self.saturationfunction(np.linspace(min(self.current), max(self.current), 1000), *params)
		
	def PlotOffsets(self, mode='linear'):
		
		if mode=='linear':
			small_data = {k:v for k,v in self.linear_params.items() if self.M.FindSiPMsize(k)=='small'}
			smalloffsets = [i[1] for i in small_data.values()]
			smallSiPMs = small_data.keys()
   
			large_data = {k:v for k,v in self.linear_params.items() if self.M.FindSiPMsize(k)=='large'}
			largeoffsets = [i[1] for i in large_data.values()]
			largeSiPMs = large_data.keys() 
   
		# elif mode=='inv_exp':
		# 	offsets = [i[-11] for i in self.saturation_params.values()]
		# 	SiPMs = self.saturation_params.keys()
	
		fig, ax = plt.subplots(figsize=(10,6))
		
		ax.set_xlabel('SiPM number')
		ax.set_ylabel('QDC offset [a.u]')
		ax.set_title(f"QDC expected for a current of 0 nA\n{self.M.titledict[self.M.mode]} illumination")

		largeSiPMs_plot = ax.bar(largeSiPMs, largeoffsets, color='orange', edgecolor='black', label='Large SiPMs')
		smallSiPMs_plot = ax.bar(smallSiPMs, smalloffsets, color='green', edgecolor='black', label='Small SiPMs')
		ax.legend()
  
		outfilename=f"{self.M.sndswpath}analysis-plots/qdc_calibration/qdcoffsets_{self.M.titledict[self.M.mode]}.png"
		fig.savefig(outfilename)
		print(f'Offsets plot save to: {outfilename}')
		
	def WriteOutOffsets(self):
		offsets_filename = f'{self.M.sndswpath}Results/Offsets/qdcoffsets_{self.M.titledict[self.M.mode]}.csv'
		# offsets_df = pd.DataFrame(self.linear_params, columns=['SiPM number', 'Linear fit'])
		# offsets_df.to_csv(offsets_filename, index=False)
		
	def OverlaySiPMs(self,bar=5):

		sipmtypekey={0:'large', 1:'small'}
		SiPMs = list(range((bar-1)*8+1, bar*8+1))

		fig, axes = plt.subplots(2,1, figsize=(10,14))
		# fig_title=fig.suptitle(f'Overlay of large and small SiPM QDC responses to laser')	
		# fig_title.set_y(0.5)
		
		for i in range(2):
			axes[i].set_title(f"Bar {bar}: QDC v current for {sipmtypekey[i]} SiPMs\n{self.M.titledict[self.M.mode]} illumination")
			axes[i].set_xlabel(self.int_title)
			axes[i].set_ylabel('QDC [a.u]')
			axes[i].grid(True, linestyle='--', alpha=0.7, color='gray')

		for idx,SiPM in enumerate(SiPMs):
			
			x=self.SetSiPM(SiPM)
			if x==-999: print(f'No runs for SiPM {SiPM}')
			# current, qdcs, errors = self.GetQDCvIntensity()

			# Plot the first scatter plot with error bars (blue dots)
			if self.M.SiPMsize=='large': axes[0].errorbar(self.current.array, self.qdcs, yerr=self.errors, color=self.colours[idx], label=f'SiPM {SiPM}')
			elif self.M.SiPMsize=='small': axes[1].errorbar(self.current.array, self.qdcs, yerr=self.errors, color=self.colours[idx], label=f'SiPM {SiPM}')

		for i in range(2): axes[i].legend()

		# Create directories if not existing already
		d=f'{self.M.sndswpath}/analysis-plots/qdc_calibration/'
		os.makedirs(d, exist_ok=True)

		plotlocation = d+f'bar{bar}_{self.M.titledict[self.M.mode]}.png'
		fig.savefig(plotlocation, bbox_inches='tight')
		print(f'Plot saved to {plotlocation}')

	def CompareRunTypes(self, style='overlay'):

		if style=='overlay':
			fig, ax = plt.subplots(figsize=(10,6))

			ax.set_title(f"QDC v {self.laser_mode} for SiPM {self.M.SiPM}, bar {self.M.bar}\nwith different illumination modes")
			ax.set_xlabel(self.int_title)
			ax.set_ylabel('QDC [a.u]')
			ax.grid(which='major',axis='both',linestyle='--')
			for idx,runtype in enumerate(self.M.runs.keys()):
				# For SiPMs with no 3-bar illumination data:
				if len(self.M.runs[runtype])==0: continue

				self.M.mode=runtype
				self.GetQDCvIntensity()

				if runtype=='SiPM':	label = f"Single SiPM {self.M.SiPM}"
				elif runtype==1:	label = f"1-bar illuminated"
				else:				label = f"{runtype}-bars illuminated"
				
				ax.errorbar(self.current.array, self.qdcs, yerr=self.errors, linestyle='-', color=self.colours[idx], label=label)
			ax.legend()	

		if style=='subplots':
			fig, axes = plt.subplots(2,2, figsize = (12, 12))
			fig_title=fig.suptitle(f"QDC v {self.laser_mode} for SiPM {self.M.SiPM}\nwith different illumination modes")

			for idx,runtype in enumerate(self.M.runs.keys()):
				self.M.mode=runtype
				# current, qdcs, errors = self.GetQDCvIntensity()
				if runtype=='SiPM':	axes[0,0].set_title(f"Single SiPM")
				elif runtype==1:	axes[0,1].set_title(f"1-bar illuminated")
				else:				axes[idx//2, idx%2].set_title(f"{runtype}-bars illuminated")

				axes[idx//2, idx%2].set_xlabel(laser_mode)
				axes[idx//2, idx%2].set_ylabel('QDC [a.u]')
				axes[idx//2, idx%2].grid(which='major', axis='both', linestyle='--')
				axes[idx//2, idx%2].errorbar(self.current.array, self.qdcs, yerr=self.errors, fmt='o', color=self.colours[idx], label=f'SiPM {self.M.SiPM}')

		# Create directories if not existing already
		d=f'{self.M.sndswpath}analysis-plots/qdc_calibration/illuminations/'
		os.makedirs(d, exist_ok=True)

		plotlocation = d+f'{style}-SiPM{self.M.SiPM}.png'
		fig.savefig(plotlocation, bbox_inches='tight')
		print(f'Plot saved to {plotlocation}')		

