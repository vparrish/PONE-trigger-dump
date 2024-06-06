#/usr/local/bin
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
parser.add_argument("-w", "--window", default=5,#ns
                    help="length of coincidence time window")
parser.add_argument("-m", "--moduleReq", default=2,
                    help="Minimum number of modules which must have PEs for an event to be considered")

args = parser.parse_args()

#not sure if this is the best way to do this
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

def findModuleMultiplicity(modules , timeWindow, req_mult):
    triggers_all_frames = []
    for frame in modules:
        triggers = []
        for mk, pmts in frame.items():
            maxMult = 0
            print("module key is "+str(mk))
            print("pmts is "+str(pmts))
            
            #this for loop gets rid of pulses with charges that are less than 0.25/too small 
            for pmt, pulses, in pmts:
                #pmt is the pmt, pulses is a pulsequeue list of I3recopulses
                #this loop gets rid of pulses that are too small- unsure if it works rn and this is something to look at later
                for p in pulses:
                    if(p.charge <0.25):
                        print("boo charge too small, removing pulse")
                        pulses.pulses.remove(p)
                    #i am not sure if this works but it doesn't actually get to this loop rn
                    if(len(pulses) == 0):
                        pmts.remove(pmt)
                        print("i have removed the pmt, all pulses too small")
            
            trigger = ModuleTrigger(mk, 0 ,0)
            startTime = np.inf
            mult = 0
                #print(pmts)
            #goes thru each pulse and determines if it's the leadtube (has the earliest starttime)
            for pmt, pulses in pmts:
                #print(pulses)
                
                for p in pulses: 
                    if p.time < startTime:
                        leadTube = pmt
                        startTime = p.time 
                    print("pulse time is "+str(p.time))
                    
                    if p.time < startTime + timeWindow:
                        mult = mult+ 1

                    if pmt==leadTube:
                        assert p.time<startTime+timeWindow
                        assert mult>0
                    #print("leadtube is "+str(leadTube))
                    #print("starttime is now "+str(startTime))
                #print("multiplicity = "+str(mult))
            #the way that assess_muon_triggers.cpp is set up rn it only reports the highest mult trigger not all triggers...
            #want to report all here so for all pmts we're just getting the multiplicities and the start times
            if mult > trigger.multiplicity:
                trigger.multiplicity = mult
                #print("multiplicity = "+str(trigger.multiplicity))
                trigger.time = startTime
                #print("multiplicity = "+str(trigger.time))
            if trigger.multiplicity >= req_mult:
                #print(trigger.multiplicity)
                triggers.append(trigger)
        print("new frame")
        triggers_all_frames.append(triggers)
    print(len(triggers))
    return triggers  
        


#we'll run this at the end and it'll return some fun arrays
#this will include the above functions 
def getData():
    frames = []
    infiles = sorted(glob(args.infile))
    for file in infiles:
        print(file)
        infile = dataio.I3File(file)
        #for i in range(0, 1):
        #    frames.append(infile.pop_daq())
        #i += 1

        while infile.more():
            #for i in range(0,1):
            try:
                frames.append(infile.pop_daq())
            except:                    
                continue
            #i += 1
    
        infile.close()
    print(len(frames))

    
    modules_div = []
    for frame in frames:
        modules = defaultdict(list)
        pulsemap = frame['PMTResponse_nonoise']
        for omkey, p in pulsemap:
            #create a dictionary with (string, om) as keys and [pmt, pulses] as items
            key = (omkey[0], omkey[1])     
            modules[key].append((omkey[2], p))
        modules_div.append(modules)
    #print(modules_div)

    #now that we have some kind of data structure for the modules, lets implement the one function
    #for i in modules:
     #   print(modules[i][0])
        #print(j)
    triggers = findModuleMultiplicity(modules_div, 10, 2)
    for i in range(len(triggers)):
        print("multiplicity= "+str(triggers[i].multiplicity))
        print("module= "+str(triggers[i].module))
        print("time= "+str(triggers[i].time))
        print("next")




getData()

        #do some more stuff that includes the functions and then returns some kind of data frame
        #not sure what to output things as tho 

#export data as some typeeeee 

