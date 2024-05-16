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



#add arguments for input, output, coincidence time window requirements and number of modules to be considered an event
parser = argparse.ArgumentParser(
    description="collects all the l0 triggers and outputs to some format not yet determined :)")

parser.add_argument("-i", "--infile", default="./test_input.i3",
                    help="input .gz files")
parser.add_argument("-o", "--outfile", default="./test_output.i3",
                    help="Write output to OUTFILE (.i3{.gz} format)")
parser.add_argument("-w", "--window", default=5,#ns
                    help="length of coincidence time window")
parser.add_argument("-m", "--moduleReq", default=2,
                    help="Minimum number of modules which must have PEs for an event to be considered")

args = parser.parse_args()

#some fun lil classes we'll need ? 
class ModuleTrigger:
    def __init__(self, module, multiplicity, time):
        self.module = module
        self.multiplicity = multiplicity
        self.time = time



#make! some! functions! 

def findModuleMultiplicity(modules , timeWindow):
    triggers = []
    for omkey,pmts, in modules:
        maxMult = 0
        for pmt, pulses, in pmts:
            if(pulses.next().GetCharge()<0.25):
                pulses.advance()
                if(pulses.empty()):
                    pmts.erase(pmt)
        #this needs to be a specific structure 
        trigger = ModuleTrigger(omkey, 0, 0)
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
    print("i have the infiles")
    for file in infiles:
        print(file)
        infile = dataio.I3File(file)
        for i in range(0,10):
            try:
                print("correct")
                print(infile.pop_DAQ())
                frames.append(infile.pop_DAQ()[i])
            except:
                #print("oops continuing")
                continue
        infile.close()
    print(frames)
    #test just with the find multiplicity module 
    for frame in frames:
        pulsemap = frame['PMTResponse_nonoise']
        modules = []
        pulseQueue = []
        for omkey, p in pulsemap:
            print("doing the thing")
            key = omkey(p[0].GetString(), p[0].GetOM())
            modules[key].update({p[0].GetPMT(): pulseQueue(p[1])})
        print(key)





getData()

        #do some more stuff that includes the functions and then returns some kind of data frame
        #not sure what to output things as tho 

#export data as some typeeeee 

