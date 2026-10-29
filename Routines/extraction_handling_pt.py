from scipy.optimize import brentq
from scipy.optimize import curve_fit
import scipy.fft as fft
import numpy as np
from datetime import datetime
import time, sys, os, pyvisa, subprocess
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.colors as mcolors
import matplotlib.cm as cm
from matplotlib.collections import LineCollection
import time
from matplotlib.colors import LinearSegmentedColormap


#Extracting collected data for analysis and appropriate plotting

sys.path.append("..")

from Routines.vna_analysis import *


def vapor_pressure_equation(temp):
    return np.power(10, (3.75641 - 95.599/(temp-1.503)))




triple_point = 24.553

pressure_triple_point = vapor_pressure_equation(triple_point)

critical_point = 44.4

pressure_critical_point = vapor_pressure_equation(critical_point)




def dataframe_loader_pressure(data_dir, low, high):



    for file in os.listdir(data_dir)[low:high]:
        read_val = (file.split('_')[-2])
        if read_val == "frequency":
            freq_df = pd.read_csv(data_dir + "\\" + file)
        elif read_val == "phaserad":
            phase_rad_df = pd.read_csv(data_dir + "\\" + file)
        elif read_val == "phase":
            phase_df = pd.read_csv(data_dir + "\\" + file)
        elif read_val == "S21":
            s21_df = pd.read_csv(data_dir + "\\" + file)
        elif read_val == "temperature":
            temp_df = pd.read_csv(data_dir + "\\" + file)
        elif read_val == "time":
            time_df = pd.read_csv(data_dir + "\\" + file)
        elif read_val == "pressure":
            pressure_df = pd.read_csv(data_dir + "\\" + file)

    freq_array = freq_df.iloc[:,1].to_numpy()
    phase_array = phase_df.iloc[:, 1:].to_numpy()
    phase_rad_array = phase_rad_df.iloc[:, 1:].to_numpy()
    temp_array = temp_df.iloc[:, -1].to_numpy()
    s21_array = s21_df.iloc[:, 1:].to_numpy()
    time_array = time_df.iloc[:, 1].to_numpy()
    pressure_array = pressure_df.iloc[:, 1].to_numpy()

    true_S_lin = 10**((s21_array - np.max(s21_array))/10)


    return freq_array, phase_array, phase_rad_array, temp_array, s21_array, true_S_lin, time_array, pressure_array



def dataframe_loader(data_dir, low, high):



    for file in os.listdir(data_dir)[low:high]:
        read_val = (file.split('_')[-2])
        if read_val == "frequency":
            freq_df = pd.read_csv(data_dir + "\\" + file)
        elif read_val == "phaserad":
            phase_rad_df = pd.read_csv(data_dir + "\\" + file)
        elif read_val == "phase":
            phase_df = pd.read_csv(data_dir + "\\" + file)
        elif read_val == "S21":
            s21_df = pd.read_csv(data_dir + "\\" + file)
        elif read_val == "temperature":
            temp_df = pd.read_csv(data_dir + "\\" + file)
        elif read_val == "time":
            time_df = pd.read_csv(data_dir + "\\" + file)

    freq_array = freq_df.iloc[:,1].to_numpy()
    phase_array = phase_df.iloc[:, 1:].to_numpy()
    phase_rad_array = phase_rad_df.iloc[:, 1:].to_numpy()
    temp_array = temp_df.iloc[:, -1].to_numpy()
    s21_array = s21_df.iloc[:, 1:].to_numpy()
    time_array = time_df.iloc[:, 1].to_numpy()

    true_S_lin = 10**((s21_array - np.max(s21_array))/10)


    return freq_array, phase_array, phase_rad_array, temp_array, s21_array, true_S_lin, time_array



def dataframe_plotter(freq_array, temp_array, s21_array, time_array, pressure_array):

    true_S_lin = 10**((s21_array - np.max(s21_array))/10)

    #fig, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4)

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3)

    a1 = ax1.plot(time_array[:], pressure_array, linewidth = 0.5, color = 'green')
    a2 = ax2.plot(time_array, temp_array, linewidth = 0.5, color = 'red')
    a3 = ax3.pcolor(time_array[:], freq_array[:], true_S_lin, cmap = 'viridis')
    #a4 = ax4.pcolor(time_array, freq_array, phase_array)

    #ax1.set_ylim([4, 20])
    ax1.set_ylabel("Pressure (mbar)", fontsize = 12)
    ax2.set_ylabel("Temperature (K)")
    ax1.set_xlabel("Time (s)", fontsize = 12)
    ax2.set_xlabel("Time (s)", fontsize = 12)
    ax3.set_xlabel("Time (s)", fontsize = 12)
    #ax4.set_xlabel("time (s)")
    #fig.colorbar(a2, label = "S$_2$$_1$ (dB)")
    fig.colorbar(a3, label = "Linear transmission")
    #fig.colorbar(a4, label = "Phase (rad)")

    fig.suptitle("Liquid Deposition", fontsize = 15)
    fig.set_size_inches([10, 4])
    fig.tight_layout()
    plt.show()


def res_frequencies(freq_array, s21_array, true_S_lin, time_array):

    f0_dep_list = []

    for slice in range (len(true_S_lin.transpose())):
        f0, q = single_scan_analysis_silent(freq_array, s21_array.transpose()[slice], true_S_lin.transpose()[slice], time_stamp = 1, forplot = False)

        f0_dep_list.append(f0)

    np.array(f0_dep_list)
    plt.plot(time_array, f0_dep_list, color = "black")
    plt.xlabel("Time (s) over deposition")
    plt.ylabel("Resonant Frequency (GHz)")

    plt.title("Frequency shift over time")

    return f0_dep_list




def baseline_loadandfit(dir2):
    fitted_f0, fitted_Q, T_data = fq_load_plot(dir2)

    temp_freq_param, temp_freq_cov = fq_fitting(T_data= T_data[:], fitted_f0 = fitted_f0[:], num_points = 501)

    return (temp_freq_param, temp_freq_cov)


def normalized_frequency_extraction(temp_array, time_array, temp_freq_param, f0_dep_list):   

    temp_adjusted_f = poly3(temp_array, *temp_freq_param)

    plt.plot(time_array, temp_adjusted_f, color = 'orange')
    plt.ticklabel_format(useOffset = False, style = 'plain')
    plt.xlabel("time (s)")
    plt.ylabel("Resonant Frequency (GHz)")
    plt.title("Temperature-induced frequency shifts")

    plt.show()

    plt.plot(time_array[:], f0_dep_list[:], color = "black")
    plt.xlabel("Time (s)", fontsize = 12)
    plt.ylabel("Resonant Frequency (GHz)", fontsize = 12)
    plt.plot(time_array[:], temp_adjusted_f[:], color = 'orange')
    plt.title("Frequency shift over time", fontsize = 15)

    plt.legend(['Measured Frequency', 'Temperature-induced frequency shift'])

    normalized_freq_shift = ((f0_dep_list - temp_adjusted_f)/temp_adjusted_f)*100

    plt.show()


    plt.plot(time_array, normalized_freq_shift, linewidth = 3)
    plt.gca().yaxis.set_major_formatter(ticker.ScalarFormatter(useOffset=False))
    plt.ticklabel_format(style = 'plain', axis = 'y')
    plt.xlabel("Time(s)", fontsize = 12)
    plt.title("Neon-induced Frequency Shift", fontsize = 15)
    plt.ylabel("Resonant Frequency shift \nas a percentage (%)")
    plt.show()

    return normalized_freq_shift



class measuredData:

    def __init__(self, filepath, elements, pres = True):

        self.filepath = filepath

        if pres:
            self.freq, self.phase, self.phase_rad, self.temp, self.s21, self.lin, self.time, self.pressure = dataframe_loader_pressure(filepath, elements[0], elements[1])
        else:
            self.freq, self.phase, self.phase_rad, self.temp, self.s21, self.lin, self.time = dataframe_loader(filepath, elements[0], elements[1])


    def plotter(self):
        dataframe_plotter(self.freq, self.temp, self.s21, self.time, self.pressure)

    def res_frequency_extraction(self):
        f0_dep_list = []

        for slice in range (len(self.lin.transpose())):
            f0, q = single_scan_analysis_silent(self.freq, self.s21.transpose()[slice], self.lin.transpose()[slice], time_stamp = 1, forplot = False)

            f0_dep_list.append(f0)

        np.array(f0_dep_list)
        plt.plot(self.time, f0_dep_list, color = "black")
        plt.xlabel("Time (s) over deposition")
        plt.ylabel("Resonant Frequency (GHz)")

        plt.title("Frequency shift over time")

        self.res_frequencies = f0_dep_list


    def freq_peak_extract(self):
        #Called if extract resonant frequencies by max value 
        res_frequencies = []
        
        for slice in range (len(self.lin.transpose())):
            res_frequencies.append(self.freq[np.argmax(self.lin.transpose()[slice])])

        self.res_frequencies = res_frequencies

    def normalized_frequency_extraction(self, baseline_params):
        temp_adjusted_f = poly3(self.temp, *baseline_params)

        plt.plot(self.time, temp_adjusted_f, color = 'orange')
        plt.ticklabel_format(useOffset = False, style = 'plain')
        plt.xlabel("time (s)")
        plt.ylabel("Resonant Frequency (GHz)")
        plt.title("Temperature-induced frequency shifts")

        plt.show()

        plt.plot(self.time[:], self.res_frequencies[:], color = "black")
        plt.xlabel("Time (s)", fontsize = 12)
        plt.ylabel("Resonant Frequency (GHz)", fontsize = 12)
        plt.plot(self.time[:], temp_adjusted_f[:], color = 'orange')
        plt.title("Frequency shift over time", fontsize = 15)

        plt.legend(['Measured Frequency', 'Temperature-induced frequency shift'])

        normalized_freq_shift = ((self.res_frequencies - temp_adjusted_f)/temp_adjusted_f)*100

        plt.show()


        plt.plot(self.time, normalized_freq_shift, linewidth = 3)
        plt.gca().yaxis.set_major_formatter(ticker.ScalarFormatter(useOffset=False))
        plt.ticklabel_format(style = 'plain', axis = 'y')
        plt.xlabel("Time(s)", fontsize = 12)
        plt.title("Neon-induced Frequency Shift", fontsize = 15)
        plt.ylabel("Resonant Frequency shift \nas a percentage (%)")
        plt.show()

        self.normalized_frequencies = normalized_freq_shift

    def triplepoint_ratefinder(self, triple_point = 24.53, rate_range = 1):

        closest_index = int(np.argmin(np.abs(self.temp - triple_point)))
        print("Closest index:", self.temp[closest_index])

        rate = abs(
        (self.temp[closest_index + rate_range] - self.temp[closest_index - rate_range])
        / (self.time[closest_index + rate_range] - self.time[closest_index - rate_range])
    )

        rate *= 60

        self.rate = rate
        return rate
    
    def min_extractor(self, filter_plot = False):
        from scipy.signal import find_peaks, savgol_filter
        smooth_temp = savgol_filter(self.temp, window_length = 11, polyorder = 3)

        minimas, _ = find_peaks(-smooth_temp, prominence = 0.1)

        min_temp = self.temp[minimas]

        if filter_plot == True:
            plt.plot(self.temp, label = 'raw')
            plt.plot(smooth_temp, label = 'smooth')
            plt.xlabel('index')
            plt.ylabel('Temp (K)')
            plt.legend()
            plt.show()

        
        self.min_frequencies = np.array(self.res_frequencies)[minimas]
        self.min_temps = min_temp
        
        return minimas, min_temp
    
    def max_extractor(self, filter_plot = False):
        from scipy.signal import find_peaks, savgol_filter
        smooth_temp = savgol_filter(self.temp, window_length = 11, polyorder = 3)

        maximas, _ = find_peaks(smooth_temp, prominence = 0.1)

        max_temp = self.temp[maximas]

        if filter_plot == True:
            plt.plot(self.temp, label = 'raw')
            plt.plot(smooth_temp, label = 'smooth')
            plt.xlabel('index')
            plt.ylabel('Temp (K)')
            plt.legend()
            plt.show()

        
        self.min_frequencies = np.array(self.res_frequencies)[maximas]
        self.max_temps = max_temp
        
        return maximas, max_temp


class measuredIQnoise:
    def __init__(self, chopping_freq, filepath):
        self.chopping_freq = chopping_freq

        self.IQvoltage = np.loadtxt(filepath, delimiter= ',', skiprows= 1)

        self.I = self.IQvoltage[:, 0]

        self.Q = self.IQvoltage[:, 1]

        self.num_points = len(self.I)

        self.time = np.linspace(0, len(self.I))/chopping_freq

    def plot_slice(self, slice):
        plt.plot(self.time[0:slice], self.I[0:slice], label = "I")
        plt.plot(self.time[0:slice],self.Q[0:slice], label = "Q")
        plt.xlabel("Time (s)")
        plt.ylabel("Voltage (V)")
        plt.legend()
        plt.show()

    def return_phase(self, unwrapped = True):
        if unwrapped:
            self.phase = np.arctan(self.Q/self.I)
        else:
            self.phase = np.arctan2(self.Q, self.I)
        print(unwrapped)

    def return_magnitude(self):
        self.mag = np.sqrt(self.I**2 + self.Q**2)

    def phase_fft(self):
        self.yf_phase = np.abs(fft.rfft(self.phase))
        self.xf_phase = fft.rfftfreq(len(self.phase), 1/self.chopping_freq)

    def magnitude_fft(self):
        self.yf_mag = np.abs(fft.rfft(self.mag))
        self.xf_mag = fft.rfftfreq(len(self.mag), 1/self.chopping_freq)

    def magnitude_psd(self):
        self.psd_mag = (self.yf_mag**2)/(len(self.mag)*self.chopping_freq)
    
    def phase_psd(self):
        self.psd_phase = (self.yf_phase**2)/(len(self.phase)*self.chopping_freq)

