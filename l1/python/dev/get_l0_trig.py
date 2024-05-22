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


#def create_pulses():
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
class PulseQueue:
    def __init__(self, pulses):
        self.it = iter(pulses)
        self.end = len(pulses)
        self.pulses = pulses

    def empty(self):
        return self.it == self.end

    def next_pulse(self):
        return next(self.it)
    
    #this function is buggy and does not work properly
    def advance(self):
        while self.it != self.end:
            #self.it += 1
            pulse =self.pulses[i]
            print(pulse)
            if i != self.end and pulse.charge >= 0.25:
                break
"""

#make! some! functions! 

def findModuleMultiplicity(modules , timeWindow, req_mult):
    triggers = []
    for mk, pmts in modules.items():
        maxMult = 0 
        #this for loop gets rid of pulses with charges that are less than 0.25/too small 
        for pmt, pulses, in pmts:
            #pmt is the pmt, pulses is a pulsequeue object
            #for each pulse in the pulsequeue object
            for p in pulses.pulses:
                
                if(p.charge <0.25):
                    print("boo charge too small, removing pulse")
                    pulses.pulses.remove(p)
                #i am not sure if this works but it doesn't actually get to this loop rn
                if(len(pulses.pulses) == 0):
                    #print("i am erasing")
                    pmts.remove(pmt)
                    print("i have removed the pmt, all pulses too small")
        
        trigger = ModuleTrigger(mk, 0 ,0)
        i = 0
        #while len(pmts) != 0:
        #for i in range(len(pmts)):
        startTime = np.inf
            #print(pmts)
        mult = 0
        for pmt, pulses in pmts:
            for p in pulses.pulses: 
                if p.time < startTime:
                    leadTube = pmt
                    startTime = p.time
                    #print("leadtube is "+str(leadTube))
                    #print("starttime is now "+str(startTime))
                    if p.time < startTime + timeWindow:
                        mult = mult+ 1

                    if pmt==leadTube:
                        assert p.time<startTime+timeWindow
                        assert mult>0
        #the way that assess_muon_triggers.cpp is set up rn it only reports the highest mult trigger not all triggers...
        #want to report all here so for all pmts we're just getting the multiplicities and the start times
        if mult > trigger.multiplicity:
            trigger.multiplicity = mult
            trigger.time = startTime
        if trigger.multiplicity:
            #print(trigger.multiplicity)
            triggers.append(trigger)

    print(len(triggers))
    #removes triggers with multiplicity < 2
    for trig in triggers:
        if trig.multiplicity < req_mult:
            triggers.remove(trig)
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
        while infile.more():
            try:
                frames.append(infile.pop_daq())
            except:
                continue
        infile.close()
    print(len(frames))
    #test just with the find multiplicity module 
    #i think we don't want to use a dictionary actually ????
    modules = defaultdict(list)
    #modules = []
    for frame in frames:
        #not sure if I should use PMTResponse with or without noise...without noise gives an assertion error
        #probably related to negative times
        pulsemap = frame['PMTResponse_nonoise']
        for omkey, p in pulsemap:
            #create a dictionary with (string, om) as keys and [pmt, pulses] as items
            key = (omkey[0], omkey[1])
            temp = PulseQueue(p)
            #so now the pulses are a PulseQueue object         
            modules[key].append((omkey[2], temp))
    #print(modules)

    #now that we have some kind of data structure for the modules, lets implement the one function
    #for i in modules:
     #   print(modules[i][0])
        #print(j)
    triggers = findModuleMultiplicity(modules, 10, 2)




getData()

        #do some more stuff that includes the functions and then returns some kind of data frame
        #not sure what to output things as tho 

#export data as some typeeeee 

