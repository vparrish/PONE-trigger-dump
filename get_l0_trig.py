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
#class ModuleTrigger:
#    def __init__(self, module, multiplicity, time):
#        self.module = module
#        self.multiplicity = multiplicity
#        self.time = time
#idk if this is right either tbh
ModuleTrigger = defaultdict(lambda:defaultdict(list))



#make! some! functions! 

def findModuleMultiplicity(modules , timeWindow):
    triggers = []
    for omkey in modules:
        maxMult = 0
        pmt = modules[omkey][0]
        print(pmt)
        pmts = modules[omkey]
        print(pmts)
        #I think I'm misunderstanding what the modules map contains in it?? 
        for pmt, pulses, in pmts:
            if(pulses.next().GetCharge()<0.25):
                pulses.advance()
                if(pulses.empty()):
                    pmts.erase(pmt)
        #this needs to be a specific structure 
        trigger = ModuleTrigger[("omkey", 0)].append(omkey)
        print(trigger)
        while not pmts.empty():
            startTime = np.inf
            for pmt, pulses in pmts:
                if pulses.next().GetTime() < startTime:
                    leadTube = pmt
                    startTime = pulses.next().GetTime()

            mult = 0
            for pmt, pulses in pmts:
                if pulses.next().GetTime() < startTime + timeWindow:
                    mult += 1
                if pmt==leadTube:
                    assert pulses.next().GetTime()<startTime+timeWindow
                    assert mult>0
            
            if mult > trigger.multiplicity:
                trigger.multiplicity = mult
                trigger.time = startTime
            pmts[leadTube].advance()
            if pmts[leadTube].empty():
                pmts.erase(leadTube)

        if trigger.multiplicity:
            triggers.append(trigger)
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
    modules = defaultdict(lambda: defaultdict(list))
    for frame in frames:
        pulsemap = frame['PMTResponse_nonoise']
        
        pulseQueue = []
        i = 0

        for omkey, p in pulsemap:
            #print(omkey)
            #print(p)
            modules[omkey] = [omkey[2], {p[0], p[-1]}]
            #print(omkey[2])
    #print(modules)

    #now that we have some kind of data structure for the modules, lets implement the one function
    #for i in modules:
     #   print(modules[i][0])
        #print(j)
    triggers = findModuleMultiplicity(modules, 10)




getData()

        #do some more stuff that includes the functions and then returns some kind of data frame
        #not sure what to output things as tho 

#export data as some typeeeee 

