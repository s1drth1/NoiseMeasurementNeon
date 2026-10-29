
"""
Configure and capture homodyne measurement data from the alazar ats9870 digital acquisition card. Refer
to Instruments.atsapi (native python wrapper) for flags/options/settings

@author: sidharth

01/18/2026

"""

from scipy.optimize import brentq
import numpy as np
from datetime import datetime
import time, sys, os, pyvisa, subprocess
import matplotlib.pyplot as plt
import time, shutil


sys.path.append("..")

from Instruments.atsapi import *
from Instruments.srs_dg535_ttlsource import *
from Instruments.SRS_SG396_lo import *
from Instruments.agilent_33220a_awg import *

global alazarboard

global bytesPerSample

global num_channels

bytesPerSample = 1

num_channels = 2

def configure_board(sample_rate, sample_rate_compute, external_clock = True, range_Q = INPUT_RANGE_PM_2_V, range_I = INPUT_RANGE_PM_2_V, trigger_external = True, trigger_timeout = 0, trigger_delay = 0):

    #configure board settings, run before every data acquisition sequence as a reset tool. Default 
    #configuration is for an external clock and trigger

    
    alazarboard = Board()
    
    if external_clock:
    
        alazarboard.setCaptureClock(EXTERNAL_CLOCK_10MHz_REF, int(1E9), CLOCK_EDGE_RISING, int(1e9/sample_rate_compute))

    else:
        alazarboard.setCaptureClock(INTERNAL_CLOCK, sample_rate, CLOCK_EDGE_RISING, 0)

    alazarboard.inputControlEx(CHANNEL_A, DC_COUPLING, range_Q, IMPEDANCE_50_OHM) #Channel 2, DC coupling, input range, 50 ohm impedance

    alazarboard.inputControlEx(CHANNEL_B, DC_COUPLING, range_I, IMPEDANCE_50_OHM)

    if trigger_external:
        alazarboard.setTriggerOperation(TRIG_ENGINE_OP_J, TRIG_ENGINE_J, TRIG_EXTERNAL, TRIGGER_SLOPE_POSITIVE, 140, TRIG_ENGINE_K, TRIG_DISABLE, TRIGGER_SLOPE_POSITIVE, 128)
        alazarboard.setExternalTrigger(DC_COUPLING, ETR_5V)

    else:
        alazarboard.setTriggerOperation(TRIG_ENGINE_OP_J, TRIG_ENGINE_J, TRIG_CHAN_A, TRIGGER_SLOPE_POSITIVE, 140, TRIG_ENGINE_K, TRIG_DISABLE, TRIGGER_SLOPE_POSITIVE, 128)
    
    
    alazarboard.setTriggerTimeOut(trigger_timeout)

    alazarboard.setTriggerDelay(trigger_delay)

    return alazarboard

def configure_trigger_mixer_settings(trig_addr = "GPIB0::15::INSTR", mixI_addr = "GPIB0::21::INSTR", mixQ_addr = "GPIB0::10::INSTR", lo_addr = "GPIB0::27::INSTR", lo_freq =2.2289E9, dc_I = 0.0132, dc_Q = -0.061, phaseQ = 49.00, chop_freq = int(1E3), chop_amplitude = 1.0, burst_cycles = 30000, channel = 4, trigger_rate = 0.03125, trigger_mode = 0, trigger_load = 0, trigger_spacing = 100E-6):

    #Configured for a square chopping wave from the Q port of the mixer

    lo_set_frequency(lo_addr, lo_freq)

    lo_set_power(lo_addr, 10.00)

    lo_on_or_off(lo_addr, "ON")

    #Q mixer for burst sequence
    awg_set_mode(mixQ_addr, "SQU")

    awg_set_amplitude(mixQ_addr, chop_amplitude)

    awg_set_frequency(mixQ_addr, chop_freq)

    awg_set_offset(mixQ_addr, dc_Q)

    awg_phaseoffset(mixQ_addr, phaseQ)

    awg_burst(mixQ_addr, "ON", burst_cycles)

    awg_trigger_config(mixQ_addr, "EXT", "POS")

    awg_on_or_off(mixQ_addr, "ON")

    #I mixer for DC offset
    awg_set_mode(mixI_addr, "DC")

    awg_set_offset(mixI_addr, dc_I)

    awg_on_or_off(mixI_addr, "ON")
    
    #Set trigger
    TTL_pulse_configuration(trig_addr, channel, trigger_rate, trigger_mode, trigger_load, trigger_spacing)

    print("Sweep settings configured")
    


def acquire_data_triggeredstreaming(alazarboard, filepath, buffers = 150, num_samples = int(10E6), acquisition_time_sec = 30, sample_rate_compute = int(100E6), to_save = True, to_plot = False):

    if to_save:
        os.makedirs(filepath, exist_ok= True)

    bytesPerBuffer = num_samples*bytesPerSample*num_channels

    samplesPerAcquisition = int(sample_rate_compute*acquisition_time_sec)

    buffersPerAcquisition = samplesPerAcquisition//num_samples

    print(samplesPerAcquisition, buffersPerAcquisition)

    bufflist = []
    for buff in range(buffers):
        bufflist.append(DMABuffer(alazarboard.handle, c_uint8, bytesPerBuffer))


    time.sleep(1)

    #enables both channes in dualport autodma traditional mode

    alazarboard.beforeAsyncRead(CHANNEL_A | CHANNEL_B,
                        0,
                        num_samples,
                        1,
                        0x7FFFFFFF,
                        ADMA_EXTERNAL_STARTCAPTURE | ADMA_TRIGGERED_STREAMING)

    for buffer in bufflist:
        alazarboard.postAsyncBuffer(buffer.addr, buffer.size_bytes) 

    start = time.perf_counter() # Keep track of when acquisition started
    try:
        alazarboard.startCapture() # Start the acquisition
        print("Capturing %d buffers. Press <enter> to abort" %
                buffersPerAcquisition)
        buffersCompleted = 0
        bytesTransferred = 0
        while (buffersCompleted < buffersPerAcquisition and not
                enter_pressed()):
            # Wait for the buffer at the head of the list of available
            # buffers to be filled by the board.
            buffer = bufflist[buffersCompleted % len(bufflist)]
            alazarboard.waitAsyncBufferComplete(buffer.addr, timeout_ms= 40000)
            buffersCompleted += 1
            bytesTransferred += buffer.size_bytes

            out = buffer.buffer
        
            if to_save:
                out.tofile(filepath + "\\"+ f"data{buffersCompleted}.bin")

            if to_plot:
                plt.plot(out)
            

            # Add the buffer to the end of the list of available buffers.
            alazarboard.postAsyncBuffer(buffer.addr, buffer.size_bytes)

            duration = time.perf_counter() - start
    finally:
        alazarboard.abortAsyncRead()


def read_and_downsample(folderpath, new_path, filename, sample_rate = int(100E6), chopping_freq = int(1E3), num_samples = int(10E6), offset = 8000, bounds_Q = 2, bounds_I = 2):

    downsampledIcomplete = []
    downsampledQcomplete = []

    
    step = int(sample_rate/chopping_freq)
    for file in os.listdir(folderpath):
        path = folderpath + '\\' + file

        data = np.fromfile(path, dtype = np.uint8).reshape([2, num_samples])

        I_data_raw = data[1, :]
        Q_data_raw = data[0, :]

        I_data = to_voltage(I_data_raw, bounds_I)
        Q_data = to_voltage(Q_data_raw, bounds_Q)

        for i in range(0, len(I_data), step):
            I_data_chunk_a = I_data[i + offset:i + int(step/2) - offset]
            I_data_chunk_b = I_data[i+int(step/2) + offset:i + int(step) - offset]

            Q_data_chunk_a = Q_data[i + offset:i + int(step/2) - offset]
            Q_data_chunk_b = Q_data[i+int(step/2) + offset:i + int(step) - offset]

            I_downsampled = return_chunk_average(I_data_chunk_a, I_data_chunk_b)
            Q_downsampled = return_chunk_average(Q_data_chunk_a, Q_data_chunk_b)

            downsampledIcomplete.append(I_downsampled)
            downsampledQcomplete.append(Q_downsampled)



    plt.plot(downsampledIcomplete)
    plt.plot(downsampledQcomplete)

    out = np.column_stack([downsampledIcomplete, downsampledQcomplete])

    os.makedirs(new_path, exist_ok= True)
    np.savetxt(new_path + "\\" + filename, out, delimiter= ',', header= 'I, Q', fmt = "%5f")

    shutil.rmtree(folderpath)

def to_voltage(data, bounds):
    return(data/255.0)*(2*bounds) - bounds

def return_chunk_average(chunk_a, chunk_b):

    chunk_a_average = np.average(chunk_a)

    chunk_b_average = np.average(chunk_b)

    return chunk_a_average - chunk_b_average

