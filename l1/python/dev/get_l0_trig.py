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

@dataclass
class ModuleTrigger:
    def __init__(self, module, multiplicity, time):
        self.module = module
        self.multiplicity = multiplicity
        self.time = time
def removal(pmts, leadTube):
    #res = [tuple(ele for ele in sub if ele != leadTube) for sub in pmts]
    for p in pmts:
        if p[0] == leadTube:
            pmts.remove(p)
    return(pmts)


#make! some! functions! 
modules = defaultdict(list)
def findModuleMultiplicity(modules , timeWindow, req_mult):

    triggers = []
    
    for mk, pmts in modules.items():
        maxMult = 0
        #print("module key is "+str(mk))
        #print("pmts is "+str(pmts))
                
        pmt_list = []      
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
                pmt_list.append(pmt)
            #print(pmt_list)   
        trigger = ModuleTrigger(mk, 0 ,0)       
       
        #print(pmts)
        #goes thru each pulse and determines if it's the leadtube (has the earliest starttime)
        
        while(len(pmts)!= 0):
            startTime = np.inf 
            
            for pmt, pulses in pmts:
                    print(len(pmts))
                    for p in pulses: 
                        print(p.time)
                        print(startTime)
                        if p.time < startTime:
                                leadTube = pmt
                                startTime = p.time 
                                #print("pulse time is "+str(p.time))
                                print("leadtube is "+str(leadTube))
                        #leadtube_list.append(leadTube)
                                print("starttime is now "+str(startTime))
            mult = 0
            for pmt, pulses in pmts:
                    #print(p ulses)
                pulse_count = 0
                for p in pulses:
                    print("length "+str(len(pulses)))
                    pulse_count +=1
                    if pulse_count < 2:
                        if p.time < startTime + timeWindow:
                            mult +=1 
                        #add this to prevent double counting for multiple pulses on one pmt? 
                        
                    print("mult is "+str(mult))
                            
                    
            if mult > trigger.multiplicity:
                trigger.multiplicity = mult
                #print("multiplicity = "+str(trigger.multiplicity))
                trigger.time = startTime
            
            #print(len(pmts))
            removal(pmts, leadTube)
            #print(len(pmts))
            
        #print("I have exited the while loop")
        if trigger.multiplicity >= req_mult:
            triggers.append(trigger)
        #print("mult is "+str(trigger.multiplicity))
        #print("mk is "+str(trigger.module))
        
    #print("new mk")        
        
    return triggers  

#function to order triggers by time
def sortTime(ModuleTrig):
    return ModuleTrig.time
def time_order(triggers):
    triggers.sort(key=sortTime)
    #print(triggers)
    return triggers
    

def getData(frame):
    modules = defaultdict(list)
    pulsemap = frame['PMTResponse_nonoise']
    for omkey, p in pulsemap:
        #create a dictionary with (string, om) as keys and [pmt, pulses] as items
        key = (omkey[0], omkey[1])     
        modules[key].append((omkey[2], p))
    #print(modules)
    triggers = findModuleMultiplicity(modules, args.window, args.moduleReq)
    #print(len(triggers))
    #for i in triggers:
        #print("mult is "+str(i.multiplicity))
        #print("time is "+str(i.time))
        #print("mk is "+str(i.module))
    
    ordered_triggers = time_order(triggers)
    #print("mult is "+str(triggers[0].multiplicity))
    #print("mk is "+str(triggers[0].module))
    #print("time is "+str(triggers[0].time))
    print("length is " +str(len(ordered_triggers)))
    for i in ordered_triggers:
        print("mk is "+str(i.module))
        print("mult is "+str(i.multiplicity))
        print("time is "+str(i.time))
        
    print("new frame")
    #return ordered_triggers



#t = I3Tray()

#t.AddModule("I3Reader", Filename= args.infile)
#t.AddModule(getData, "getData", Streams = [icetray.I3Frame.DAQ])

#t.Execute()


