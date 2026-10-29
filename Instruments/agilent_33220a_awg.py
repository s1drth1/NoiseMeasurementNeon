"""
Code to configure a awg source, this is specifically written for the AG 33220a

Created: 01/09/2026

@author: sidharthduthaluru

"""

import pyvisa, time

import numpy

rm = pyvisa.ResourceManager()


def awg_set_mode(address, mode):
    #Set function mode (DC, square, sine, etc.)
    awg = rm.open_resource(address)

    awg.write(f"SOUR:FUNC {mode}")

    time.sleep(0.1)

    print(awg.query("SOUR:FUNC?"))

    time.sleep(0.1)

    awg.close()


def awg_set_frequency(address, freq):
    #Setting frequency
    awg = rm.open_resource(address)

    awg.write(f"SOUR:FREQ {freq}")

    time.sleep(0.1)

    print(awg.query("SOUR:FREQ?"))
    
    time.sleep(0.1)

    awg.close()

def awg_set_amplitude(address, amplitude_peaktopeak):
    #Set signal amplitude
    awg = rm.open_resource(address)

    awg.write(f"SOUR:VOLT {amplitude_peaktopeak}")

    time.sleep(0.1)

    print(awg.query("SOUR:VOLT?"))

    time.sleep(0.1)

    awg.close()

def awg_set_offset(address, offset):
    #Set signal DC offset
    awg = rm.open_resource(address)

    awg.write(f"SOUR:VOLT:OFFS {offset}")

    time.sleep(0.1)

    print(awg.query("SOUR:VOLT:OFFS?"))

    time.sleep(0.1)

    awg.close()

def awg_on_or_off(address, on_or_off):
    #Toggle output
    awg = rm.open_resource(address)

    awg.write(f"OUTP {on_or_off}")

    time.sleep(0.1)

    print(awg.query("OUTP?"))

    time.sleep(0.1)

    awg.close()

def awg_phaseoffset(address, phase):
    #Configure phase offset, important for clock synchronization
    awg = rm.open_resource(address)

    awg.write(f"PHAS {phase}")
    time.sleep(0.1)
    print(awg.query("PHASE?"))
    time.sleep(0.1)
    awg.close()

def awg_burst(address, state, ncycles):
    #Configure a triggered signal burst
    awg = rm.open_resource(address)
    awg.write(f"BURST:STAT {state}")
    awg.write(f"BURST:NCYCLES {ncycles}")
    time.sleep(0.1)
    print(awg.query("BURST:STAT?"))
    print(awg.query("BURST:NCYCLES?"))
    time.sleep(0.1)
    awg.close()
    
def awg_trigger_config(address, source, slope):
    #Configure trigger used before burst signal 
    awg = rm.open_resource(address)

    awg.write(f"TRIG:SOUR {source}")
    awg.write(f"TRIG:SLOP {slope}")
    time.sleep(0.1)
    print(awg.query("TRIG:SOUR?"))
    print(awg.query("TRIG:SLOP?"))
    time.sleep(0.1)
    awg.close()