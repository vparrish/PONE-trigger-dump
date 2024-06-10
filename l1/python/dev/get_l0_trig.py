#!/usr/bin/env python
'''
date: 6 March 2024
author: jmgarriz

Script to retrieve L0 events from simulation. Loose structure based on get_data_tables.py
script from Jean-Pierre

'''
from I3Tray import *
from icecube import icetray, dataio, dataclasses
from icecube import phys_services
from icecube.icetray import I3Units

from glob import glob
import pandas as pd
import argparse
import os
import numpy as np
from collections import defaultdict
from dataclasses import dataclass



#add arguments for input, output, coincidence time window requirements and number of modules to be considered an event
parser = argparse.ArgumentParser(
    description="collects all the l0 triggers and outputs to some format not yet determined :)")

parser.add_argument("-i", "--infile", default="/mnt/research/IceCube/PONE/jp_pone_sim/pmtsim/GenerateSingleMuons_39_pmtsim.i3.zst",
                    help="input .gz files")
parser.add_argument("-o", "--outfile", default="./out.i3",
                    help="Write output to OUTFILE (.i3{.gz} format)")
parser.add_argument("-w", "--window", default=10,#ns
                    help="length of coincidence time window")
parser.add_argument("-m", "--moduleReq", default=2,
                    help="Minimum number of modules which must have PEs for an event to be considered")

args = parser.parse_args()
#t = I3Tray()
#t.context['I3RandomService'] = phys_services.I3GSLRandomService(42)
#t.Add('I3Reader', FilenameList=args.infile)
#

@dataclass
class ModuleTrigger:
    def __init__(self, module, multiplicity, time):
        self.module = module
        self.multiplicity = multiplicity
        self.time = time

"""
this pulseQueue class doesn't actually need to be here
class PulseQueue:
    def __init__(self, pulses):
        self.pulses = pulses
        self.index = 0
    
    def empty(self):
        return self.index >= len(self.pulses)
    
    def next_pulse(self):
        #if self.empty():
        #    raise StopIteration("No more pulses")
        return self.pulses[self.index]
    
    def advance(self):
        while self.index < len(self.pulses):
            self.index += 1
            if self.index != len(self.pulses) and self.pulses[self.index].get_charge() >= 0.25:
                break        
"""

#make! some! functions! 
modules = defaultdict(list)
def findModuleMultiplicity(modules , timeWindow, req_mult):
    #triggers_all_frames = []
    #for frame in modules:
    triggers = []
    for mk, pmts in modules.items():
        maxMult = 0
        #print("module key is "+str(mk))
        #print("pmts is "+str(pmts))
                
                #this for loop gets rid of pulses with charges that are less than 0.25/too small 
        for pmt, pulses, in pmts:
                #pmt is the pmt, pulses is a pulsequeue list of I3recopulses
                #this loop gets rid of pulses that are too small- unsure if it works rn and this is something to look at later
            for p in pulses:
                if(p.charge <0.25):
                    print("boo charge too small, removing pulse")
                    pulses.remove(p)
                    #i am not sure if this works but it doesn't actually get to this loop rn
                if(len(pulses) == 0):
                    pmts.remove(pmt)
                    print("i have removed the pmt, all pulses too small")
                
        trigger = ModuleTrigger(mk, 0 ,0)       
        mult = 0
                    #print(pmts)
                #goes thru each pulse and determines if it's the leadtube (has the earliest starttime)
                #maybe instead of this loop we can just order by time so that then we sequentially got thru and let each pmt
                #be the lead tube on its own 
            #leadtube_list = []
            #leadTube = 999
        startTime = np.inf 
        for pmt, pulses in pmts:
            for p in pulses: 
                if p.time < startTime:
                        leadTube = pmt
                        startTime = p.time 
                        #print("pulse time is "+str(p.time))
                #print("leadtube is "+str(leadTube))
                #print("starttime is now "+str(startTime))

        for pmt, pulses in pmts:
                #print(pulses)
            for p in pulses:
                if p.time < startTime + timeWindow:
                    #print(startTime)
                    mult = mult+ 1
        if mult > trigger.multiplicity:
            trigger.multiplicity = mult
            #print("multiplicity = "+str(trigger.multiplicity))
            trigger.time = startTime

        if trigger.multiplicity >= req_mult:
                #print(trigger.multiplicity)
            triggers.append(trigger)
        #triggers_all_frames.append(triggers)
        #print("new frame")
            
        #print(len(triggers_all_frames))
    return triggers  

"""
#i don;t know if this works
def time_order(triggers):
    for frame in triggers:
        for i in range(1, len(frame)):
            key_item = frame.time[i]
            print(key_item)
            j = i - 1
            while j >= 0 and frame.time[j] > key_item:
            # Shift the value one position to the left
            # and reposition j to point to the next element
            # (from right to left)
                frame.time[j + 1] = frame.time[j]
                j -= 1

        # When you finish shifting the elements, you can position
        # `key_item` in its correct location
            frame.time[j + 1] = key_item

    return triggers
"""

def getData(frame):
    modules = defaultdict(list)
    pulsemap = frame['PMTResponse_nonoise']
    for omkey, p in pulsemap:
        #create a dictionary with (string, om) as keys and [pmt, pulses] as items
        key = (omkey[0], omkey[1])     
        modules[key].append((omkey[2], p))
    print(modules)
    triggers = findModuleMultiplicity(modules, args.window, args.moduleReq)
    #print("mult is "+str(triggers[0].multiplicity))
    #print("mk is "+str(triggers[0].module))
    #print("time is "+str(triggers[0].time))
    print("length is " +str(len(triggers)))
    for i in triggers:
        print("mult is "+str(i.multiplicity))
        print("time is "+str(i.time))
        print("mk is "+str(i.module))
#
t = I3Tray()
i = 0
for i in range(0,1):
    t.AddModule("I3Reader", Filename= args.infile)
    t.AddModule(getData, "getData", Streams = [icetray.I3Frame.DAQ])

    t.Execute()
    i +=1


