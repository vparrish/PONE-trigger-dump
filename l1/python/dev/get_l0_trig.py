#!/usr/bin/env python
'''
date: 6 March 2024
author: jmgarriz

Script to retrieve L0 events from simulation. 
!!WORK IN PROGRESS!! 
'''
from I3Tray import *
from icecube import icetray, dataio, dataclasses
from icecube import phys_services
from icecube.icetray import I3Units

from glob import glob
import argparse
import os
import numpy as np
from collections import defaultdict
from dataclasses import dataclass
import l1_lc_alg as l1



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

