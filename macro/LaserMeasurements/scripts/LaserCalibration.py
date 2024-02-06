import csv,os,json
import numpy as np 
import pandas as pd
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
plt.rcParams["text.usetex"]
# plt.rcParams['xtick.labelsize'] = 16
# plt.rcParams['ytick.labelsize'] = 16
# plt.rcParams['axes.labelsize'] = 16
# plt.rcParams['figure.titlesize'] = 18
plt.rcParams['font.size'] = 16
import matplotlib.patches as patches

class LaserCalibration(object):
    def __init__(self):
        self.path=f'/afs/cern.ch/user/a/aconsnd/LaserMeasurements'
        self.datafilename=f"/afs/cern.ch/user/a/aconsnd/LaserMeasurements/configuration/avg_adc_extrapolation_external_rising_edge_10kHz_run1run2run3_ch30.csv"
        self.GetData()

        self.config=self.ReadConfig()
        self.intensityprofile=self.GetIntensityProfile()
        self.fig, self.axes = plt.subplots(3, 1, figsize=(10,20)) # plot: adc v uncalibrated, adc v calibrated, n_gamma v calibrated 
        self.calibratedADC_2_incidentPhotons_per_mmsq = 1/self.config['calibration_channel_gain'] * 1/self.config['calibration_channel_area'] * 1/self.config['calibration_channel_PDE']
        self.linearstartpoint=20

    def ReadConfig(self):
        configfile = f"{self.path}/configuration/config.json"
        with open(configfile, 'r') as f:
            d = json.load(f) 
        return d        

    def GetIntensityProfile(self):
        with open(f'{self.path}/configuration/intensityprofile.csv', 'r') as f:
            d={i[0]:i[1:] for i in csv.reader(f)}
        for i in d: d[i] = [float(j) for j in d[i]]
        return d        

    def GetData(self):
        df = pd.read_csv(self.datafilename, index_col=False, header=0)
        # df = df.dropna(subset=['Average adc [a.u]'])
        df = df[df['Average adc [a.u]']>0]
        self.df=df

    def PlotIntensityvCurrent(self):
        
        # self.axes[0].errorbar(np.array(self.df["Laser intensity [%]"]), np.array(self.df["Average adc [a.u]"]), yerr=np.array(self.df['Error [a.u]']), color='C0', marker='o', label='Channel 30')
        self.axes[0].errorbar(self.df["Laser intensity [%]"].array, self.df["Average adc [a.u]"].array, yerr=self.df['Error [a.u]'].array, color='C0', marker='o', label='Channel 30')

        # Get line of best fit
        self.FitLine()
        # red_chi2 = self.GetRedChi2()
        self.axes[0].errorbar(self.x_fit, self.y_fit, label='Linear fit', color='orange', linestyle='-')

        self.axes[0].set_xlabel('Laser intensity [%]')
        self.axes[0].set_ylabel('Average adc [a.u]')
        self.axes[0].grid(True)
        self.axes[0].set_title('Uncalibrated intensities')
        self.axes[0].legend()

        # Draw fit on plot
        # textbox = patches.Rectangle((12,62), 25,5, linewidth=1, edgecolor='black', facecolor='white')
        # plt.gca().add_patch(textbox)
        # eqn_text = f'y={self.coefficients[0]:.2f}x {self.coefficients[1]:.2f}'
        # chi2_text = f'$\\chi^2_\\nu$ = {red_chi2:.2f}' 
        # self.axes[0].text(55,500, eqn_text, fontsize=12, color='black')
        # self.axes[0].text(55,250, chi2_text, fontsize=12, color='black')
        
    def linearfunction(self, x, a, b):	
        return a*x + b

    def FitLine(self):

        # I want to do a straight line fit between 20% and 100%, then just linear interpolation between 0% and 20%
        linear_range = int(np.where(self.df["Laser intensity [%]"]==self.linearstartpoint)[0])
        linear_x, linear_y, adc_errors = self.df["Laser intensity [%]"][linear_range:].array, self.df["Average adc [a.u]"][linear_range:].array, self.df["Error [a.u]"][linear_range:].array

        m_0 = (linear_y[2] - linear_y[1]) / (linear_x[2] - linear_x[1])
        c_0 = linear_y[2] - linear_x[2]*m_0
        self.linear_params, self.linear_covariance = curve_fit(self.linearfunction, linear_x, linear_y, p0=[m_0, c_0], sigma=adc_errors)

        self.x_fit = np.linspace(linear_x[0], linear_x[-1], 1000)
        self.y_fit = self.x_fit*self.linear_params[0] + self.linear_params[1]
    
    def GetRedChi2(self):

        residuals = self.df["Average adc [a.u]"] - np.polyval(self.linear_params, self.df["Laser intensity [%]"])
        chi_2 = np.sum(residuals**2)
        degree=1
        red_chi2 = chi_2/(degree+1)
        return red_chi2

    def GetCalibratedIntensities(self):

        # self.calibrated_currents={i:0 for i in self.x}
        intensities = sorted(list(set([i for s in self.intensityprofile.values() for i in s])))
        intensities = list(set(intensities + list(self.df['Laser intensity [%]'])))

        self.calibrated_data = pd.DataFrame(columns=['Laser intensity [%]','Calibrated laser intensity [%]','Average adc [a.u]','Average adc error [a.u]','Incident photons / mm_sq'])

        # Interpolate data less than 20% intensity
        interpolation_range = self.df.loc[self.df["Laser intensity [%]"] == self.linearstartpoint]
        row_numbers = interpolation_range.index.to_list()   # Get row number where intensity = 20%
        interp_df = self.df.iloc[:row_numbers[0]+1]     # Copy dataframe up to 20%

        linear_data = self.df.iloc[row_numbers[0]:]     # Linear fit range is 20% and above

        # print(f'Len intensities: {len(intensities)}')

        for i in range(len(intensities)):
            intensity = intensities[i]

            if intensity<=self.linearstartpoint:

                adc_interpolated = np.interp(intensity, interp_df["Laser intensity [%]"], interp_df["Average adc [a.u]"])
                
                # Use interpolated current and fitted straight line to determine the expected intensity
                intensity_exp = 1/self.linear_params[0]*(adc_interpolated-self.linear_params[1])
                
                # Get n_incident photons per mm_sq by diving by calibration channel gain, PDE and surface area. 
                n_incident_photons = adc_interpolated*self.calibratedADC_2_incidentPhotons_per_mmsq

                # Get error on ADC either from data file or estimate using the relative uncertainty of the next data point
                adc_relerror = self.GetADCerror(intensity)
                adc_error = adc_interpolated*adc_relerror

                new_row=[intensity,intensity_exp,adc_interpolated,adc_error,n_incident_photons] 
                self.calibrated_data.loc[i] = new_row

            else:   
                adc_interpolated = np.interp(intensity, linear_data["Laser intensity [%]"], linear_data["Average adc [a.u]"])
                n_incident_photons = adc_interpolated*self.calibratedADC_2_incidentPhotons_per_mmsq

                adc_relerror = self.GetADCerror(intensity)
                adc_error = adc_interpolated*adc_relerror

                new_row=[intensity,intensity,adc_interpolated,adc_error,n_incident_photons] 
                # new_row=[intensity,intensity,adc_interpolated,n_incident_photons] 
                self.calibrated_data.loc[i] = new_row                

        # sorted_indices = np.argsort(self.calibrated_data[:, 0])
        # self.calibrated_data[sorted_indices]
        self.calibrated_data = self.calibrated_data.sort_values(by='Laser intensity [%]')

    def GetADCerror(self, intensity):
        
        # print(intensity)
        if intensity in self.df['Laser intensity [%]'].array: 
            adc_relerror = self.df.loc[self.df['Laser intensity [%]'] == intensity, 'Error [a.u]'].values[0] / self.df.loc[self.df['Laser intensity [%]'] == intensity, 'Average adc [a.u]'].values[0]


        else:
            # index of next largest value of laser intensity 
            next_largest_index = self.df[self.df['Laser intensity [%]'] > intensity].idxmin()['Laser intensity [%]']
            print(f'Next largest index: {next_largest_index}')
            # relative uncertainty on adc of next largest point
            adc_relerror = self.df.loc[next_largest_index, 'Error [a.u]'] / self.df.loc[next_largest_index, 'Average adc [a.u]']

        return adc_relerror

    def PlotCalibratedIntensities(self):

        # self.axes[1].errorbar(self.calibrated_data[:, 1], self.calibrated_data[:, 2], yerr=self.df['Error [a.u]'].array, marker='o', label='Calibrated intensity')
        self.axes[1].errorbar(self.calibrated_data['Calibrated laser intensity [%]'].array, self.calibrated_data['Average adc [a.u]'].array, yerr=self.calibrated_data['Average adc error [a.u]'].array, color='C0', marker='o', label='Calibrated intensity\nchannel 30')
        x_fit = np.linspace(min(self.calibrated_data['Calibrated laser intensity [%]'].array),max(self.calibrated_data['Calibrated laser intensity [%]'].array),1000)
        y_fit = self.linear_params[0]*x_fit + self.linear_params[1]
        self.axes[1].plot(x_fit,y_fit, label='Linear fit', color='C1', linestyle='-')

        self.axes[1].set_xlabel('Calibrated laser intensity [%]')
        self.axes[1].set_ylabel('Average adc [a.u]')
        self.axes[1].grid(True)
        self.axes[1].set_title('Calibrated intensities')
        self.axes[1].legend() 
        
    def PlotNincidentPhotons(self):

        # self.axes[2].errorbar(self.calibrated_data[:, 1], self.calibrated_data[:, 3], yerr=np.sqrt(self.calibrated_data[:, 3]), color='C1', label='Incident $n_{\gamma}$ per $mm^{2}$')
        
        ax2 = self.axes[2].twinx()
        
        # detectedphotons = self.calibrated_data['Incident photons / mm_sq'] / self.cal
        # self.axes[2].errorbar(self.calibrated_data['Calibrated laser intensity [%]'], detectedphotons, yerr=np.sqrt(detectedphotons), color='C0', label='Detected $n_{\gamma}$')
        # self.axes[2].plot(self.x_fit, self.y_fit, label='Linear fit', color='orange', linestyle='-')

        self.axes[2].set_xlabel('Calibrated laser intensity [%]')
        self.axes[2].set_ylabel('$n_{\gamma} / mm^{2}$')
        ax2.set_ylabel('$n_{\gamma}$ on full sensor')
        self.axes[2].grid(True)
        self.axes[2].set_title('Number of incident photons on calibration channel')

        self.axes[2].errorbar(self.calibrated_data['Calibrated laser intensity [%]'].array, self.calibrated_data['Incident photons / mm_sq'].array, yerr=np.sqrt(self.calibrated_data['Incident photons / mm_sq']), color='C1', label='Incident $n_{\gamma}$ per $mm^{2}$')
        
        # Plot incident photons on full SiPM
        ph_SiPM, counting_error = 0.74*36*self.calibrated_data['Incident photons / mm_sq'].array, 0.74*36*np.sqrt(self.calibrated_data['Incident photons / mm_sq'])
        ax2.errorbar(self.calibrated_data['Calibrated laser intensity [%]'].array, ph_SiPM, yerr=counting_error, color='C1', label='Incident $n_{\gamma}$ on large US SiPM')
        # self.axes[2].legend()            

        plotloc=f'/afs/cern.ch/user/a/aconsnd/LaserMeasurements/configuration/lasercalibration.png'
        plt.savefig(plotloc,bbox_inches='tight')
        print(f'Plot save to {plotloc}')

    def WriteOutCalibratedIntensities(self):
        calibrated_intensities_file = f"{self.path}/configuration/calibrated_intensities.csv"
        self.calibrated_data.to_csv(calibrated_intensities_file)   
        # with open(calibrated_intensities_file, 'w', newline='') as f:
        #     writer=csv.writer(f)
        #     writer.writerow(["Laser intensity [%]", "Calibrated laser intensity [%]", "Average adc [a.u]", "Incident photons / mm_sq"])

        #     for x in self.calibrated_data:
        #         writer.writerow(x)

        print(f'Calibrated intensities written to {calibrated_intensities_file}')

LC=LaserCalibration()
LC.PlotIntensityvCurrent()
LC.GetCalibratedIntensities()
LC.PlotCalibratedIntensities()
LC.PlotNincidentPhotons()
LC.WriteOutCalibratedIntensities()
        
