"""
Configuration/control file for the SRS DG535 TTL source to be used as a trigger for
synchronized measurement

Note that these controls don't follow the standard SCPI commandset. Therefore this file is highly
specific towards this particular usecase

Refer to the SRS DG535 programming manual for further modifications
@author: sidharth

Date: 01/17/2026
"""


import pyvisa, time

import numpy as np



rm = pyvisa.ResourceManager()

def TTL_pulse_configuration(address, channel = 4, trigger_rate = 0.03125, trigger_mode = 0, trigger_load = 0, trigger_spacing = 100E-6):
    #Specifically written for output on A/B --- need to update further

    ttl = rm.open_resource(address)

    ttl.write(f"TM {trigger_mode}; TR 0,{trigger_rate}; TZ {channel},{trigger_load}; OM {channel},0; DT 2,1,{trigger_spacing}; DT 3,2,{trigger_spacing}")
    print(f"TM {trigger_mode}; TR 0,{trigger_rate}; TZ {channel},{trigger_load}; OM {channel},0; DT 2,1,{trigger_spacing}; DT 3,2,{trigger_spacing}")
    print(f"{trigger_spacing}")
    

    time.sleep(0.1)

    print("TTL confirmed")

    ttl.close()